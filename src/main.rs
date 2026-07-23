use anyhow::{Context, Result, bail};
use clap::Parser;
use dicom_core::Tag;
use dicom_object::{FileDicomObject, InMemDicomObject, open_file};
use dicom_pixeldata::PixelDecoder;
use num_traits::NumCast;
use rayon::prelude::*;
use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use walkdir::WalkDir;

mod nifti_write;
mod spatial;

use spatial::SliceInfo;

#[derive(Debug, Parser)]
#[command(version)]
struct Cli {
    /// Input directory containing DICOM files
    input: PathBuf,

    /// Output directory for NIfTI files
    #[arg(short, long)]
    output: PathBuf,

    /// Number of parallel threads (0 = auto)
    #[arg(short = 'j', long, default_value = "0")]
    threads: usize,
}

#[derive(Debug, Clone)]
struct DicomSliceInfo {
    path: PathBuf,
    relative_dir: PathBuf,
    series_uid: String,
    series_number: i32,
    series_description: String,
    study_date: String,
    rows: u16,
    cols: u16,
    bits_allocated: u16,
    pixel_representation: u16,
    slice_info: SliceInfo,
    subvol_key: i64,
    subvol_label: String,
    acquisition_number: i32,
    echo_number: i32,
    instance_number: i32,
    image_type: String,
    is_derived: bool,
    is_phase: bool,
    is_real: bool,
    is_imaginary: bool,
    is_magnitude: bool,
    content_time: f64,
    trigger_delay_time: f64,
    rows_cols_key: u64,
    // TCGA-style metadata
    patient_id: String,
    patient_name: String,
    study_instance_uid: String,
    study_description: String,
    modality: String,
    manufacturer: String,
    manufacturer_model: String,
    magnetic_field_strength: f64,
    repetition_time: f64,
    echo_time: f64,
    flip_angle: f64,
    protocol_name: String,
    body_part: String,
    institution_name: String,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct GroupKey {
    relative_dir: PathBuf,
    series_uid: String,
    subvol_key: i64,
}

fn read_dicom_header(path: &Path, input_root: &Path) -> Result<DicomSliceInfo> {
    let obj = open_file(path).with_context(|| format!("Opening {}", path.display()))?;

    let series_uid = get_string(&obj, Tag(0x0020, 0x000E)).unwrap_or_default();
    let series_number = get_i32(&obj, Tag(0x0020, 0x0011)).unwrap_or(0);
    let series_description = get_string(&obj, Tag(0x0008, 0x103E)).unwrap_or_default();
    let study_date = get_string(&obj, Tag(0x0008, 0x0020)).unwrap_or_default();
    let rows = get_u16(&obj, Tag(0x0028, 0x0010)).unwrap_or(0);
    let cols = get_u16(&obj, Tag(0x0028, 0x0011)).unwrap_or(0);
    let bits_allocated = get_u16(&obj, Tag(0x0028, 0x0100)).unwrap_or(16);
    let pixel_representation = get_u16(&obj, Tag(0x0028, 0x0103)).unwrap_or(0);

    let ipp = get_f64_vec(&obj, Tag(0x0020, 0x0032)).unwrap_or_default();
    let iop = get_f64_vec(&obj, Tag(0x0020, 0x0037)).unwrap_or_default();
    let ps = get_f64_vec(&obj, Tag(0x0028, 0x0030)).unwrap_or_default();
    let st = get_f64(&obj, Tag(0x0018, 0x0050)).unwrap_or(1.0);

    let position = if ipp.len() >= 3 {
        [ipp[0], ipp[1], ipp[2]]
    } else {
        [0.0, 0.0, 0.0]
    };
    let orientation = if iop.len() >= 6 {
        [iop[0], iop[1], iop[2], iop[3], iop[4], iop[5]]
    } else {
        [1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
    };
    let pixel_spacing = if ps.len() >= 2 { [ps[0], ps[1]] } else { [1.0, 1.0] };

    let slice_info = SliceInfo {
        position,
        orientation,
        pixel_spacing,
        slice_thickness: st,
    };

    // EchoNumbers (0x0018, 0x0086)
    let echo_number = get_i32(&obj, Tag(0x0018, 0x0086)).unwrap_or(0);

    let (subvol_key, subvol_label) = determine_subvolume(&obj);
    let acquisition_number = get_i32(&obj, Tag(0x0020, 0x0012)).unwrap_or(0);
    let instance_number = get_i32(&obj, Tag(0x0020, 0x0013)).unwrap_or(0);

    // ImageType (0008,0008) — full parsing for derived/phase/real/imaginary/magnitude detection
    let image_type_raw = get_string(&obj, Tag(0x0008, 0x0008)).unwrap_or_default();
    let image_type_upper = image_type_raw.to_uppercase();
    let is_derived = image_type_upper.contains("DERIVED");
    let is_phase = image_type_upper.contains("_P_")
        || image_type_upper.contains("\\P\\")
        || image_type_upper.contains("PHASE");
    let is_real = image_type_upper.contains("_R_")
        || image_type_upper.contains("\\R\\")
        || image_type_upper.contains("_REAL_");
    let is_imaginary = image_type_upper.contains("_I_")
        || image_type_upper.contains("\\I\\")
        || image_type_upper.contains("_IMAGINARY_");
    let is_magnitude = image_type_upper.contains("_M_")
        || image_type_upper.contains("\\M\\")
        || image_type_upper.contains("_MAGNITUDE_");
    // Extract 3rd component for mDIXON (W/IP/OP/F)
    let image_type = {
        let parts: Vec<&str> = image_type_raw.split('\\').collect();
        if parts.len() >= 3 {
            parts[2].trim().to_string()
        } else {
            String::new()
        }
    };

    // ContentTime (0008,0033) — for time-based phase grouping
    let content_time = get_string(&obj, Tag(0x0008, 0x0033))
        .and_then(|s| s.trim().parse::<f64>().ok())
        .unwrap_or(0.0);

    // TriggerDelayTime (0020,9153) or TriggerTime (0018,1060)
    let trigger_delay_time = get_f64(&obj, Tag(0x0020, 0x9153))
        .or_else(|| get_f64(&obj, Tag(0x0018, 0x1060)))
        .unwrap_or(0.0);

    // Matrix size key for dimension consistency check
    let rows_cols_key = ((rows as u64) << 16) | (cols as u64);

    // TCGA-style metadata
    let patient_id = get_string(&obj, Tag(0x0010, 0x0020)).unwrap_or_default();
    let patient_name = get_string(&obj, Tag(0x0010, 0x0010)).unwrap_or_default();
    let study_instance_uid = get_string(&obj, Tag(0x0020, 0x000D)).unwrap_or_default();
    let study_description = get_string(&obj, Tag(0x0008, 0x1030)).unwrap_or_default();
    let modality = get_string(&obj, Tag(0x0008, 0x0060)).unwrap_or_default();
    let manufacturer = get_string(&obj, Tag(0x0008, 0x0070)).unwrap_or_default();
    let manufacturer_model = get_string(&obj, Tag(0x0008, 0x1090)).unwrap_or_default();
    let magnetic_field_strength = get_f64(&obj, Tag(0x0018, 0x0087)).unwrap_or(0.0);
    let repetition_time = get_f64(&obj, Tag(0x0018, 0x0080)).unwrap_or(0.0);
    let echo_time = get_f64(&obj, Tag(0x0018, 0x0081)).unwrap_or(0.0);
    let flip_angle = get_f64(&obj, Tag(0x0018, 0x1314)).unwrap_or(0.0);
    let protocol_name = get_string(&obj, Tag(0x0018, 0x1030)).unwrap_or_default();
    let body_part = get_string(&obj, Tag(0x0018, 0x0015)).unwrap_or_default();
    let institution_name = get_string(&obj, Tag(0x0008, 0x0080)).unwrap_or_default();

    let parent = path.parent().unwrap_or(path);
    let rel = parent.strip_prefix(input_root).unwrap_or(parent);
    let relative_dir = rel.parent().unwrap_or(rel).to_path_buf();

    Ok(DicomSliceInfo {
        path: path.to_path_buf(),
        relative_dir,
        series_uid,
        series_number,
        series_description,
        study_date,
        rows,
        cols,
        bits_allocated,
        pixel_representation,
        slice_info,
        subvol_key,
        subvol_label,
        acquisition_number,
        echo_number,
        instance_number,
        image_type,
        is_derived,
        is_phase,
        is_real,
        is_imaginary,
        is_magnitude,
        content_time,
        trigger_delay_time,
        rows_cols_key,
        patient_id,
        patient_name,
        study_instance_uid,
        study_description,
        modality,
        manufacturer,
        manufacturer_model,
        magnetic_field_strength,
        repetition_time,
        echo_time,
        flip_angle,
        protocol_name,
        body_part,
        institution_name,
    })
}

fn determine_subvolume(obj: &FileDicomObject<InMemDicomObject>) -> (i64, String) {
    // B-value detection MUST come first: Philips DWI files have TemporalPositionIdentifier=1
    // for ALL slices, which would mask the b-value differentiation if checked first.

    // Strategy 1: B-value (highest priority — DWI is the most important split)
    // Only use b-value > 0 as a definitive split key; b=0 is ambiguous (present in
    // non-DWI series like DIXON dynamics). b=0 groups get labeled in post-processing.
    // Standard DiffusionBValue (0018,9087) — used by Philips, Siemens, GE
    if let Some(bv) = get_f64(obj, Tag(0x0018, 0x9087)) {
        let bv_int = bv.round() as i64;
        if bv_int > 0 {
            return (bv_int, format!("b{bv_int}"));
        }
    }
    // Also try reading as integer (some scanners store IS instead of FD)
    if let Some(bv) = get_i32(obj, Tag(0x0018, 0x9087)) {
        let bv_int = bv as i64;
        if bv_int > 0 {
            return (bv_int, format!("b{bv_int}"));
        }
    }
    // GE private (0043,1039) — only if value > 0
    if let Some(vals) = get_f64_vec(obj, Tag(0x0043, 0x1039)) {
        if !vals.is_empty() {
            let mut bv = vals[0] as i64;
            if bv > 1_000_000_000 {
                bv -= 1_000_000_000;
            }
            if bv > 0 {
                return (bv, format!("b{bv}"));
            }
        }
    }
    // Siemens private (0019,100c) — only if value > 0
    if let Some(bv) = get_f64(obj, Tag(0x0019, 0x100C)) {
        let bv_int = bv.round() as i64;
        if bv_int > 0 {
            return (bv_int, format!("b{bv_int}"));
        }
    }

    // Strategy 2: TemporalPositionIdentifier (0x0020, 0x0100) — for DCE/dynamic series
    if let Some(tp) = get_i32(obj, Tag(0x0020, 0x0100)) {
        if tp > 0 {
            return (tp as i64, format!("phase{tp}"));
        }
    }

    (0, String::new())
}

// ---- DICOM value extraction helpers ----

fn get_string(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<String> {
    obj.element(tag)
        .ok()
        .and_then(|e| e.to_str().ok().map(|s| s.trim().to_string()))
}

fn get_i32(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<i32> {
    obj.element(tag).ok().and_then(|e| e.to_int::<i32>().ok())
}

fn get_u16(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<u16> {
    obj.element(tag).ok().and_then(|e| e.to_int::<u16>().ok())
}

fn get_f64(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<f64> {
    obj.element(tag).ok().and_then(|e| {
        e.to_float64()
            .ok()
            .or_else(|| e.to_int::<i32>().ok().map(|v| v as f64))
    })
}

fn get_f64_vec(obj: &FileDicomObject<InMemDicomObject>, tag: Tag) -> Option<Vec<f64>> {
    obj.element(tag).ok().and_then(|e| {
        e.to_multi_float64().ok().or_else(|| {
            e.to_str().ok().map(|s| {
                s.split('\\')
                    .filter_map(|v| v.trim().parse::<f64>().ok())
                    .collect()
            })
        })
    })
}

fn collect_dicom_files(input: &Path) -> Vec<PathBuf> {
    WalkDir::new(input)
        .into_iter()
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_file())
        .filter(|e| {
            let p = e.path();
            let ext = p.extension().and_then(|s| s.to_str()).unwrap_or("");
            ext.eq_ignore_ascii_case("dcm")
                || ext.eq_ignore_ascii_case("dicom")
                || ext.is_empty()
                || (!ext.contains("nii")
                    && !ext.contains("json")
                    && !ext.contains("csv")
                    && !ext.contains("txt")
                    && !ext.contains("dat"))
        })
        .map(|e| e.path().to_path_buf())
        .collect()
}

#[derive(Debug)]
struct ConversionRecord {
    output_path: String,
    relative_dir: String,
    patient_id: String,
    patient_name: String,
    study_date: String,
    study_instance_uid: String,
    study_description: String,
    series_instance_uid: String,
    series_number: i32,
    series_description: String,
    modality: String,
    manufacturer: String,
    manufacturer_model: String,
    institution_name: String,
    magnetic_field_strength: f64,
    body_part: String,
    protocol_name: String,
    repetition_time: f64,
    echo_time: f64,
    flip_angle: f64,
    subvolume_label: String,
    n_dicom_files: usize,
    shape_x: usize,
    shape_y: usize,
    shape_z: usize,
    pixel_spacing_x: f64,
    pixel_spacing_y: f64,
    slice_spacing: f64,
    dicom_rows: u16,
    dicom_cols: u16,
    dicom_slice_thickness: f64,
    position_range_mm: f64,
    spacing_uniformity: f64,
    bits_allocated: u16,
}

fn convert_group(slices: &mut [DicomSliceInfo], output_dir: &Path) -> Result<ConversionRecord> {
    if slices.is_empty() {
        bail!("Empty slice group");
    }

    let ref_info = &slices[0].slice_info;
    let normal = ref_info.slice_normal();

    slices.sort_by(|a, b| {
        let da = a.slice_info.dot_normal(&normal);
        let db = b.slice_info.dot_normal(&normal);
        da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
    });

    let n_slices = slices.len();
    let rows = slices[0].rows as usize;
    let cols = slices[0].cols as usize;
    let bits = slices[0].bits_allocated;
    let signed = slices[0].pixel_representation == 1;

    let desc = &slices[0].series_description;
    let sn = slices[0].series_number;
    let date = &slices[0].study_date;
    let subvol = &slices[0].subvol_label;
    let base_name = if subvol.is_empty() {
        format!("{}_{:04}_{}", date, sn, sanitize_filename(desc))
    } else {
        format!(
            "{}_{:04}_{}_{}",
            date,
            sn,
            sanitize_filename(desc),
            subvol
        )
    };
    let patient_dir = output_dir.join(&slices[0].relative_dir);
    std::fs::create_dir_all(&patient_dir)?;
    let out_path = patient_dir.join(format!("{base_name}.nii.gz"));

    let affine = spatial::compute_affine(
        &slices[0].slice_info,
        &slices[n_slices - 1].slice_info,
        rows,
        cols,
        n_slices,
    );

    let slice_spacing =
        (affine[0][2].powi(2) + affine[1][2].powi(2) + affine[2][2].powi(2)).sqrt();

    let position_range_mm = if n_slices > 1 {
        let first_pos = slices[0].slice_info.dot_normal(&normal);
        let last_pos = slices[n_slices - 1].slice_info.dot_normal(&normal);
        (last_pos - first_pos).abs()
    } else {
        0.0
    };

    let spacing_uniformity = if n_slices > 2 {
        let gaps: Vec<f64> = (1..n_slices)
            .map(|i| {
                (slices[i].slice_info.dot_normal(&normal)
                    - slices[i - 1].slice_info.dot_normal(&normal))
                .abs()
            })
            .collect();
        let mean = gaps.iter().sum::<f64>() / gaps.len() as f64;
        if mean > 1e-10 {
            let variance =
                gaps.iter().map(|g| (g - mean).powi(2)).sum::<f64>() / gaps.len() as f64;
            variance.sqrt() / mean
        } else {
            0.0
        }
    } else {
        0.0
    };

    let volume = if bits <= 8 && !signed {
        decode_and_stack_u8(slices, rows, cols, n_slices)?
    } else if bits == 16 && !signed {
        decode_and_stack::<u16>(slices, rows, cols, n_slices)?
    } else if bits == 16 && signed {
        decode_and_stack::<i16>(slices, rows, cols, n_slices)?
    } else {
        decode_and_stack::<i16>(slices, rows, cols, n_slices)?
    };

    nifti_write::write_nifti_i16(&out_path, &volume, &affine, [rows, cols, n_slices])?;

    let s = &slices[0];
    Ok(ConversionRecord {
        output_path: out_path.to_string_lossy().to_string(),
        relative_dir: s.relative_dir.to_string_lossy().to_string(),
        patient_id: s.patient_id.clone(),
        patient_name: s.patient_name.clone(),
        study_date: s.study_date.clone(),
        study_instance_uid: s.study_instance_uid.clone(),
        study_description: s.study_description.clone(),
        series_instance_uid: s.series_uid.clone(),
        series_number: s.series_number,
        series_description: s.series_description.clone(),
        modality: s.modality.clone(),
        manufacturer: s.manufacturer.clone(),
        manufacturer_model: s.manufacturer_model.clone(),
        institution_name: s.institution_name.clone(),
        magnetic_field_strength: s.magnetic_field_strength,
        body_part: s.body_part.clone(),
        protocol_name: s.protocol_name.clone(),
        repetition_time: s.repetition_time,
        echo_time: s.echo_time,
        flip_angle: s.flip_angle,
        subvolume_label: s.subvol_label.clone(),
        n_dicom_files: n_slices,
        shape_x: cols,
        shape_y: rows,
        shape_z: n_slices,
        pixel_spacing_x: s.slice_info.pixel_spacing[1],
        pixel_spacing_y: s.slice_info.pixel_spacing[0],
        slice_spacing,
        dicom_rows: s.rows,
        dicom_cols: s.cols,
        dicom_slice_thickness: s.slice_info.slice_thickness,
        position_range_mm,
        spacing_uniformity,
        bits_allocated: s.bits_allocated,
    })
}

fn decode_and_stack<T: NumCast + Copy + Default + Send + Sync + 'static>(
    slices: &[DicomSliceInfo],
    rows: usize,
    cols: usize,
    n_slices: usize,
) -> Result<Vec<i16>> {
    let mut volume = vec![0i16; rows * cols * n_slices];

    for (slice_idx, slice_info) in slices.iter().enumerate() {
        let obj = open_file(&slice_info.path)?;
        let pixel_data = obj.decode_pixel_data()?;
        let offset = slice_idx * rows * cols;

        match pixel_data.to_vec_frame::<T>(0) {
            Ok(frame_data) => {
                for (i, &val) in frame_data.iter().enumerate() {
                    if i < rows * cols {
                        volume[offset + i] = NumCast::from(val).unwrap_or(0);
                    }
                }
            }
            Err(_) => {
                // Raw-byte fallback for LUT failures (e.g. 12-bit data with negative rescale)
                let raw_bytes = pixel_data.data();
                let n_pixels = rows * cols;
                for i in 0..n_pixels.min(raw_bytes.len() / 2) {
                    let idx = i * 2;
                    let val = u16::from_le_bytes([raw_bytes[idx], raw_bytes[idx + 1]]);
                    volume[offset + i] = val as i16;
                }
            }
        }
    }

    Ok(volume)
}

fn decode_and_stack_u8(
    slices: &[DicomSliceInfo],
    rows: usize,
    cols: usize,
    n_slices: usize,
) -> Result<Vec<i16>> {
    let mut volume = vec![0i16; rows * cols * n_slices];

    for (slice_idx, slice_info) in slices.iter().enumerate() {
        let obj = open_file(&slice_info.path)?;
        let pixel_data = obj.decode_pixel_data()?;
        let frame_data = pixel_data.to_vec_frame::<u8>(0)?;
        let offset = slice_idx * rows * cols;
        for (i, &val) in frame_data.iter().enumerate() {
            if i < rows * cols {
                volume[offset + i] = val as i16;
            }
        }
    }

    Ok(volume)
}

fn sanitize_filename(s: &str) -> String {
    s.chars()
        .map(|c| {
            if c.is_alphanumeric() || c == '-' || c == '_' || c == '+' {
                c
            } else {
                '_'
            }
        })
        .collect::<String>()
        .trim_matches('_')
        .to_string()
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.threads > 0 {
        rayon::ThreadPoolBuilder::new()
            .num_threads(cli.threads)
            .build_global()
            .ok();
    }

    eprintln!("Scanning DICOM files in {} ...", cli.input.display());
    let files = collect_dicom_files(&cli.input);
    eprintln!("Found {} candidate files", files.len());

    // Phase 1: Read all headers in parallel
    eprintln!("Reading DICOM headers ...");
    let errors = Mutex::new(Vec::new());
    let input_root = &cli.input;
    let headers: Vec<DicomSliceInfo> = files
        .par_iter()
        .filter_map(|p| match read_dicom_header(p, input_root) {
            Ok(h) => {
                if h.rows > 0 && h.cols > 0 && !h.series_uid.is_empty() {
                    Some(h)
                } else {
                    None
                }
            }
            Err(_e) => {
                errors.lock().unwrap().push(p.clone());
                None
            }
        })
        .collect();

    let n_errors = errors.lock().unwrap().len();
    eprintln!(
        "Parsed {} valid DICOM slices ({} skipped)",
        headers.len(),
        n_errors
    );

    // Phase 1.5: Pre-filter non-imaging slices
    // Remove LOCALIZER images (dcm2niix: isLocalizer detection)
    let headers: Vec<DicomSliceInfo> = headers
        .into_iter()
        .filter(|h| {
            let it = h.image_type.to_uppercase();
            // Keep everything that is not a LOCALIZER
            it != "LOCALIZER"
        })
        .collect();

    // Phase 2: Group by series + subvolume
    let mut groups: HashMap<GroupKey, Vec<DicomSliceInfo>> = HashMap::new();
    for h in headers {
        let key = GroupKey {
            relative_dir: h.relative_dir.clone(),
            series_uid: h.series_uid.clone(),
            subvol_key: h.subvol_key,
        };
        groups.entry(key).or_default().push(h);
    }

    // Collect per-series info for post-processing
    let mut series_subvols: HashMap<(PathBuf, String), std::collections::HashSet<i64>> =
        HashMap::new();
    for key in groups.keys() {
        series_subvols
            .entry((key.relative_dir.clone(), key.series_uid.clone()))
            .or_default()
            .insert(key.subvol_key);
    }

    let mut groups_iter_labels: HashMap<(PathBuf, String), Vec<String>> = HashMap::new();
    for (key, slices) in &groups {
        if let Some(s) = slices.first() {
            groups_iter_labels
                .entry((key.relative_dir.clone(), key.series_uid.clone()))
                .or_default()
                .push(s.subvol_label.clone());
        }
    }

    let mut final_groups: HashMap<String, Vec<DicomSliceInfo>> = HashMap::new();
    for (key, mut slices) in groups {
        let has_multiple_subvols = series_subvols
            .get(&(key.relative_dir.clone(), key.series_uid.clone()))
            .map(|s| s.len() > 1)
            .unwrap_or(false);

        if !has_multiple_subvols {
            for s in &mut slices {
                s.subvol_label.clear();
                s.subvol_key = 0;
            }
        } else if key.subvol_key == 0 {
            // For DWI series: label the b=0 group only if siblings are b-value based
            let is_bval_series = groups_iter_labels
                .get(&(key.relative_dir.clone(), key.series_uid.clone()))
                .map(|labels| labels.iter().any(|l| l.starts_with('b') && l != "b0"))
                .unwrap_or(false);
            if is_bval_series {
                for s in &mut slices {
                    if s.subvol_label.is_empty() {
                        s.subvol_label = "b0".to_string();
                    }
                }
            }
        }

        let group_name = if has_multiple_subvols {
            format!("{}_{}", key.series_uid, key.subvol_key)
        } else {
            key.series_uid.clone()
        };

        // For groups with no subvolume differentiation, try splitting strategies
        // Ranked by coverage frequency (most common split reasons first)
        if !has_multiple_subvols && slices.len() > 1 {
            // === STRATEGY 1: Dimension consistency (dcm2niix: isDimensionVaries) ===
            // Slices with different matrix sizes MUST NOT be stacked together.
            // This catches mixed-resolution series before any other logic.
            let distinct_dims: std::collections::HashSet<u64> =
                slices.iter().map(|s| s.rows_cols_key).collect();
            if distinct_dims.len() > 1 {
                let mut sub_groups: HashMap<u64, Vec<DicomSliceInfo>> = HashMap::new();
                for s in slices {
                    sub_groups.entry(s.rows_cols_key).or_default().push(s);
                }
                for (dim_idx, (_dim, mut sub_slices)) in sub_groups.into_iter().enumerate() {
                    let label = format!("dim{dim_idx}");
                    for s in &mut sub_slices {
                        s.subvol_label = label.clone();
                    }
                    let name = format!("{}_dim{}", group_name, dim_idx);
                    final_groups.insert(name, sub_slices);
                }
                continue;
            }

            // === STRATEGY 2: Derived vs Original (dcm2niix: isDerived) ===
            // Do not stack DERIVED and ORIGINAL images together.
            // Catches ADC, TRACEW, FA maps mixed with source DWI.
            let has_derived = slices.iter().any(|s| s.is_derived);
            let has_original = slices.iter().any(|s| !s.is_derived);
            if has_derived && has_original {
                let mut sub_groups: HashMap<bool, Vec<DicomSliceInfo>> = HashMap::new();
                for s in slices {
                    sub_groups.entry(s.is_derived).or_default().push(s);
                }
                for (is_der, mut sub_slices) in sub_groups {
                    let label = if is_der {
                        "derived".to_string()
                    } else {
                        "original".to_string()
                    };
                    for s in &mut sub_slices {
                        s.subvol_label = label.clone();
                    }
                    let name = format!("{}_{}", group_name, label);
                    final_groups.insert(name, sub_slices);
                }
                continue;
            }

            // === STRATEGY 3: Phase/Real/Imaginary/Magnitude splitting ===
            // (dcm2niix: isHasPhase/isHasReal/isHasImaginary; dicom2nifti: ImageType checks)
            // Complex-valued series produce multiple image types in same series.
            let has_multi_complex = {
                let n_types = [
                    slices.iter().any(|s| s.is_phase),
                    slices.iter().any(|s| s.is_real),
                    slices.iter().any(|s| s.is_imaginary),
                    slices.iter().any(|s| s.is_magnitude),
                ]
                .iter()
                .filter(|&&v| v)
                .count();
                n_types > 1
            };
            if has_multi_complex {
                let mut sub_groups: HashMap<String, Vec<DicomSliceInfo>> = HashMap::new();
                for s in slices {
                    let label = if s.is_phase {
                        "ph".to_string()
                    } else if s.is_real {
                        "real".to_string()
                    } else if s.is_imaginary {
                        "imaginary".to_string()
                    } else {
                        "mag".to_string()
                    };
                    sub_groups.entry(label).or_default().push(s);
                }
                for (label, mut sub_slices) in sub_groups {
                    for s in &mut sub_slices {
                        s.subvol_label = label.clone();
                    }
                    let name = format!("{}_{}", group_name, label);
                    final_groups.insert(name, sub_slices);
                }
                continue;
            }

            // === STRATEGY 4: Multi-echo splitting (T2*, DIXON, multi-echo GRE) ===
            // (dcm2niix: echoNum/TE varies → separate; dicom2nifti: not handled)
            let distinct_echos: std::collections::HashSet<i32> =
                slices.iter().map(|s| s.echo_number).collect();
            if distinct_echos.len() > 1 {
                let mut sub_groups: HashMap<i32, Vec<DicomSliceInfo>> = HashMap::new();
                for s in slices {
                    sub_groups.entry(s.echo_number).or_default().push(s);
                }
                for (echo, mut sub_slices) in sub_groups {
                    for s in &mut sub_slices {
                        s.subvol_label = format!("echo{echo}");
                    }
                    let name = format!("{}_echo{}", group_name, echo);
                    final_groups.insert(name, sub_slices);
                }
                continue;
            }

            // === STRATEGY 5: ImageType splitting (mDIXON-All: W/IP/OP/F) ===
            // (dcm2niix: handles via imageType text; dicom2nifti: filters LOCALIZER by ImageType)
            let distinct_img_types: std::collections::HashSet<&str> = slices
                .iter()
                .filter(|s| !s.image_type.is_empty())
                .map(|s| s.image_type.as_str())
                .collect();
            if distinct_img_types.len() > 1 {
                let mut sub_groups: HashMap<String, Vec<DicomSliceInfo>> = HashMap::new();
                for s in slices {
                    let key = if s.image_type.is_empty() {
                        "unknown".to_string()
                    } else {
                        s.image_type.clone()
                    };
                    sub_groups.entry(key).or_default().push(s);
                }
                for (img_type, mut sub_slices) in sub_groups {
                    let label = img_type.to_lowercase();
                    for s in &mut sub_slices {
                        s.subvol_label = label.clone();
                    }
                    let name = format!("{}_{}", group_name, label);
                    final_groups.insert(name, sub_slices);
                }
                continue;
            }

            // === STRATEGY 6: AcquisitionNumber splitting ===
            // (dcm2niix: stacks by default but #mySegmentByAcq compile flag exists)
            // (dicom2nifti Siemens: old code used AcqNum for classic 4D grouping)
            let distinct_acq: std::collections::HashSet<i32> =
                slices.iter().map(|s| s.acquisition_number).collect();
            if distinct_acq.len() > 1 && distinct_acq.len() < slices.len() {
                let mut sub_groups: HashMap<i32, Vec<DicomSliceInfo>> = HashMap::new();
                for s in slices {
                    sub_groups.entry(s.acquisition_number).or_default().push(s);
                }
                for (acq, mut sub_slices) in sub_groups {
                    for s in &mut sub_slices {
                        s.subvol_label = format!("acq{acq}");
                    }
                    let name = format!("{}_acq{}", group_name, acq);
                    final_groups.insert(name, sub_slices);
                }
                continue;
            }

            // === STRATEGY 7: Same-position stacking (CEST, localizers) ===
            // When all slices are at essentially the same spatial position but
            // have distinct AcquisitionNumbers — split by AcqNum.
            if distinct_acq.len() > 1 && slices.len() > 2 {
                let normal = slices[0].slice_info.slice_normal();
                let positions: Vec<f64> = slices
                    .iter()
                    .map(|s| s.slice_info.dot_normal(&normal))
                    .collect();
                let pos_min = positions.iter().cloned().fold(f64::INFINITY, f64::min);
                let pos_max = positions.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
                let pos_range = (pos_max - pos_min).abs();
                let slice_thick = slices[0].slice_info.slice_thickness;
                if pos_range < slice_thick.max(1.0) {
                    let mut sub_groups: HashMap<i32, Vec<DicomSliceInfo>> = HashMap::new();
                    for s in slices {
                        sub_groups.entry(s.acquisition_number).or_default().push(s);
                    }
                    for (acq, mut sub_slices) in sub_groups {
                        for s in &mut sub_slices {
                            s.subvol_label = format!("acq{acq}");
                        }
                        let name = format!("{}_acq{}", group_name, acq);
                        final_groups.insert(name, sub_slices);
                    }
                    continue;
                }
            }

            // === STRATEGY 8: TriggerDelayTime splitting (Philips cardiac) ===
            // (dcm2niix: triggerDelayTime varies → not stacked for Philips non-ASL)
            let distinct_triggers: std::collections::HashSet<i64> = slices
                .iter()
                .map(|s| (s.trigger_delay_time * 100.0).round() as i64)
                .collect();
            if distinct_triggers.len() > 1 && distinct_triggers.len() < slices.len() {
                let mut sub_groups: HashMap<i64, Vec<DicomSliceInfo>> = HashMap::new();
                for s in slices {
                    let key = (s.trigger_delay_time * 100.0).round() as i64;
                    sub_groups.entry(key).or_default().push(s);
                }
                for (trig_idx, (_trig, mut sub_slices)) in
                    sub_groups.into_iter().enumerate()
                {
                    let label = format!("trig{trig_idx}");
                    for s in &mut sub_slices {
                        s.subvol_label = label.clone();
                    }
                    let name = format!("{}_trig{}", group_name, trig_idx);
                    final_groups.insert(name, sub_slices);
                }
                continue;
            }

            // === STRATEGY 9: Repeated-position phase splitting (UIH DCE) ===
            // When no TemporalPositionIdentifier or AcquisitionNumber is available,
            // detect repeated slice locations and split by InstanceNumber-based phases.
            // (dicom2nifti generic/siemens: uses position-direction reversal detection)
            if slices.len() > 3 {
                let normal = slices[0].slice_info.slice_normal();
                let quantize = |v: f64| -> i64 { (v * 10.0).round() as i64 };
                let unique_positions: std::collections::HashSet<i64> = slices
                    .iter()
                    .map(|s| quantize(s.slice_info.dot_normal(&normal)))
                    .collect();
                let n_unique = unique_positions.len();
                let n_total = slices.len();
                if n_unique > 1 && n_total > n_unique && n_total % n_unique == 0 {
                    let n_phases = n_total / n_unique;
                    if n_phases > 1 && n_phases <= 200 {
                        let mut pos_counts: HashMap<i64, usize> = HashMap::new();
                        for s in &slices {
                            *pos_counts
                                .entry(quantize(s.slice_info.dot_normal(&normal)))
                                .or_default() += 1;
                        }
                        let all_same_count =
                            pos_counts.values().all(|&c| c == n_phases);

                        if all_same_count {
                            slices.sort_by_key(|s| s.instance_number);
                            for phase in 0..n_phases {
                                let start = phase * n_unique;
                                let end = start + n_unique;
                                let mut sub: Vec<DicomSliceInfo> =
                                    slices[start..end].to_vec();
                                let label = format!("phase{}", phase + 1);
                                for s in &mut sub {
                                    s.subvol_label = label.clone();
                                }
                                let name =
                                    format!("{}_phase{}", group_name, phase + 1);
                                final_groups.insert(name, sub);
                            }
                            continue;
                        }
                    }
                }

                // === STRATEGY 10: Direction-reversal detection ===
                // (dicom2nifti: _classic_get_grouped_dicoms / get_grouped_dicoms)
                // Sort by InstanceNumber, compute inter-slice direction vectors,
                // and split when direction reverses (indicates new volume).
                slices.sort_by_key(|s| s.instance_number);
                let mut split_points: Vec<usize> = Vec::new();
                if slices.len() >= 3 {
                    let get_pos =
                        |s: &DicomSliceInfo| -> [f64; 3] { s.slice_info.position };
                    let mut prev_dir: Option<[f64; 3]> = None;
                    for i in 1..slices.len() {
                        let cur_pos = get_pos(&slices[i]);
                        let prev_pos = get_pos(&slices[i - 1]);
                        let diff = [
                            cur_pos[0] - prev_pos[0],
                            cur_pos[1] - prev_pos[1],
                            cur_pos[2] - prev_pos[2],
                        ];
                        let norm = (diff[0].powi(2) + diff[1].powi(2) + diff[2].powi(2))
                            .sqrt();
                        if norm < 1e-6 {
                            continue;
                        }
                        let dir = [diff[0] / norm, diff[1] / norm, diff[2] / norm];
                        if let Some(pd) = prev_dir {
                            let dot =
                                pd[0] * dir[0] + pd[1] * dir[1] + pd[2] * dir[2];
                            if dot < 0.95 {
                                split_points.push(i);
                                prev_dir = None;
                                continue;
                            }
                        }
                        prev_dir = Some(dir);
                    }
                }
                if split_points.len() >= 1 {
                    // Verify all resulting groups have the same size
                    let mut boundaries: Vec<usize> = Vec::new();
                    boundaries.push(0);
                    boundaries.extend_from_slice(&split_points);
                    boundaries.push(slices.len());
                    let group_sizes: Vec<usize> = boundaries
                        .windows(2)
                        .map(|w| w[1] - w[0])
                        .collect();
                    let all_equal = group_sizes.iter().all(|&s| s == group_sizes[0]);
                    if all_equal && group_sizes[0] >= 2 && group_sizes.len() > 1 {
                        for (phase, window) in
                            boundaries.windows(2).enumerate()
                        {
                            let mut sub: Vec<DicomSliceInfo> =
                                slices[window[0]..window[1]].to_vec();
                            let label = format!("phase{}", phase + 1);
                            for s in &mut sub {
                                s.subvol_label = label.clone();
                            }
                            let name =
                                format!("{}_phase{}", group_name, phase + 1);
                            final_groups.insert(name, sub);
                        }
                        continue;
                    }
                }
            }
        }

        final_groups.insert(group_name, slices);
    }

    // Post-process step 1: split groups with different orientations (multi-plane localizers)
    let mut orientation_split_groups: HashMap<String, Vec<DicomSliceInfo>> = HashMap::new();
    for (name, slices) in final_groups {
        if slices.len() < 2 {
            orientation_split_groups.insert(name, slices);
            continue;
        }
        let quantize_ori = |s: &DicomSliceInfo| -> [i32; 6] {
            let o = &s.slice_info.orientation;
            [
                (o[0] * 100.0).round() as i32,
                (o[1] * 100.0).round() as i32,
                (o[2] * 100.0).round() as i32,
                (o[3] * 100.0).round() as i32,
                (o[4] * 100.0).round() as i32,
                (o[5] * 100.0).round() as i32,
            ]
        };
        let first_ori = quantize_ori(&slices[0]);
        let has_multi_ori = slices.iter().any(|s| quantize_ori(s) != first_ori);
        if has_multi_ori {
            let mut ori_groups: HashMap<[i32; 6], Vec<DicomSliceInfo>> = HashMap::new();
            for s in slices {
                let key = quantize_ori(&s);
                ori_groups.entry(key).or_default().push(s);
            }
            for (idx, (_ori, mut sub_slices)) in ori_groups.into_iter().enumerate() {
                for s in &mut sub_slices {
                    s.subvol_label = format!("plane{idx}");
                }
                let sub_name = format!("{}_plane{}", name, idx);
                orientation_split_groups.insert(sub_name, sub_slices);
            }
        } else {
            orientation_split_groups.insert(name, slices);
        }
    }
    let final_groups = orientation_split_groups;

    // Post-process step 2: split groups with inconsistent slice spacing
    let mut spacing_split_groups: HashMap<String, Vec<DicomSliceInfo>> = HashMap::new();
    for (name, mut slices) in final_groups {
        if slices.len() < 4 {
            spacing_split_groups.insert(name, slices);
            continue;
        }
        let normal = slices[0].slice_info.slice_normal();
        slices.sort_by(|a, b| {
            let da = a.slice_info.dot_normal(&normal);
            let db = b.slice_info.dot_normal(&normal);
            da.partial_cmp(&db).unwrap_or(std::cmp::Ordering::Equal)
        });

        let positions: Vec<f64> = slices
            .iter()
            .map(|s| s.slice_info.dot_normal(&normal))
            .collect();
        let mut gaps: Vec<f64> = Vec::new();
        for i in 1..positions.len() {
            gaps.push((positions[i] - positions[i - 1]).abs());
        }
        if gaps.is_empty() {
            spacing_split_groups.insert(name, slices);
            continue;
        }

        let mut sorted_gaps = gaps.clone();
        sorted_gaps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let median_gap = sorted_gaps[sorted_gaps.len() / 2];

        let threshold = median_gap * 5.0;
        let mut split_indices: Vec<usize> = Vec::new();
        for (i, g) in gaps.iter().enumerate() {
            if *g > threshold && median_gap > 0.01 {
                split_indices.push(i + 1);
            }
        }

        if !split_indices.is_empty() {
            let mut valid = true;
            let mut prev = 0;
            for &idx in &split_indices {
                if idx - prev < 3 {
                    valid = false;
                    break;
                }
                prev = idx;
            }
            if slices.len() - prev < 3 {
                valid = false;
            }
            if !valid {
                split_indices.clear();
            }
        }

        if split_indices.is_empty() {
            spacing_split_groups.insert(name, slices);
        } else {
            let mut start = 0;
            let mut part_idx = 0;
            let mut all_indices: Vec<usize> = split_indices;
            all_indices.push(slices.len());
            for end in all_indices {
                let mut sub: Vec<DicomSliceInfo> = slices[start..end].to_vec();
                let existing_label = sub[0].subvol_label.clone();
                let new_label = if existing_label.is_empty() {
                    format!("part{part_idx}")
                } else {
                    format!("{existing_label}_part{part_idx}")
                };
                for s in &mut sub {
                    s.subvol_label = new_label.clone();
                }
                let sub_name = format!("{}_part{}", name, part_idx);
                spacing_split_groups.insert(sub_name, sub);
                start = end;
                part_idx += 1;
            }
        }
    }
    let final_groups = spacing_split_groups;

    eprintln!("Grouped into {} series/sub-volumes", final_groups.len());

    // Phase 3: Convert each group to NIfTI
    std::fs::create_dir_all(&cli.output)?;

    let results: Vec<(String, Result<ConversionRecord>)> = final_groups
        .into_par_iter()
        .map(|(name, mut slices)| {
            let result = convert_group(&mut slices, &cli.output);
            (name, result)
        })
        .collect();

    let mut success = 0;
    let mut failed = 0;
    let mut records: Vec<ConversionRecord> = Vec::new();
    let mut error_records: Vec<(String, String)> = Vec::new();

    for (name, result) in results {
        match result {
            Ok(rec) => {
                success += 1;
                eprintln!("  OK: {}", rec.output_path);
                records.push(rec);
            }
            Err(e) => {
                failed += 1;
                eprintln!("  FAIL [{}]: {}", name, e);
                error_records.push((name, format!("{e}")));
            }
        }
    }

    let csv_path = cli.output.join("conversion_mapping.csv");
    write_mapping_csv(&csv_path, &records)?;
    eprintln!(
        "Mapping CSV → {}  ({} entries)",
        csv_path.display(),
        records.len()
    );

    if !error_records.is_empty() {
        let err_path = cli.output.join("error_log.csv");
        write_error_csv(&err_path, &error_records)?;
        eprintln!(
            "Error log → {}  ({} errors)",
            err_path.display(),
            error_records.len()
        );
    }

    eprintln!("\nDone: {success} converted, {failed} failed");
    Ok(())
}

fn write_mapping_csv(path: &Path, records: &[ConversionRecord]) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);

    writeln!(
        f,
        "output_path,relative_dir,patient_id,patient_name,study_date,\
study_instance_uid,study_description,series_instance_uid,series_number,\
series_description,modality,manufacturer,manufacturer_model,institution_name,\
magnetic_field_strength,body_part,protocol_name,repetition_time,echo_time,\
flip_angle,subvolume_label,n_dicom_files,shape_x,shape_y,shape_z,\
pixel_spacing_x,pixel_spacing_y,slice_spacing,\
dicom_rows,dicom_cols,dicom_slice_thickness,position_range_mm,spacing_uniformity,bits_allocated"
    )?;

    for r in records {
        writeln!(
            f,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{:.2},{},{},{:.2},{:.2},{:.2},{},{},{},{},{},{:.4},{:.4},{:.4},\
{},{},{:.4},{:.4},{:.6},{}",
            csv_escape(&r.output_path),
            csv_escape(&r.relative_dir),
            csv_escape(&r.patient_id),
            csv_escape(&r.patient_name),
            &r.study_date,
            csv_escape(&r.study_instance_uid),
            csv_escape(&r.study_description),
            csv_escape(&r.series_instance_uid),
            r.series_number,
            csv_escape(&r.series_description),
            &r.modality,
            csv_escape(&r.manufacturer),
            csv_escape(&r.manufacturer_model),
            csv_escape(&r.institution_name),
            r.magnetic_field_strength,
            csv_escape(&r.body_part),
            csv_escape(&r.protocol_name),
            r.repetition_time,
            r.echo_time,
            r.flip_angle,
            csv_escape(&r.subvolume_label),
            r.n_dicom_files,
            r.shape_x,
            r.shape_y,
            r.shape_z,
            r.pixel_spacing_x,
            r.pixel_spacing_y,
            r.slice_spacing,
            r.dicom_rows,
            r.dicom_cols,
            r.dicom_slice_thickness,
            r.position_range_mm,
            r.spacing_uniformity,
            r.bits_allocated,
        )?;
    }
    Ok(())
}

fn write_error_csv(path: &Path, errors: &[(String, String)]) -> Result<()> {
    use std::io::Write;
    let mut f = std::io::BufWriter::new(std::fs::File::create(path)?);
    writeln!(f, "group_key,error")?;
    for (key, err) in errors {
        writeln!(f, "{},{}", csv_escape(key), csv_escape(err))?;
    }
    Ok(())
}

fn csv_escape(s: &str) -> String {
    if s.contains(',') || s.contains('"') || s.contains('\n') {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}
