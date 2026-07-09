use anyhow::Result;
use flate2::write::GzEncoder;
use flate2::Compression;
use std::io::Write;
use std::path::Path;

/// Write a 3D i16 volume as a gzipped NIfTI-1 file.
pub fn write_nifti_i16(
    path: &Path,
    data: &[i16],
    affine: &[[f64; 4]; 4],
    dims: [usize; 3], // [rows, cols, slices]
) -> Result<()> {
    let rows = dims[0] as i16;
    let cols = dims[1] as i16;
    let slices = dims[2] as i16;

    // Pixel dimensions from affine columns
    let pixdim_x =
        (affine[0][0].powi(2) + affine[1][0].powi(2) + affine[2][0].powi(2)).sqrt() as f32;
    let pixdim_y =
        (affine[0][1].powi(2) + affine[1][1].powi(2) + affine[2][1].powi(2)).sqrt() as f32;
    let pixdim_z =
        (affine[0][2].powi(2) + affine[1][2].powi(2) + affine[2][2].powi(2)).sqrt() as f32;

    let mut header = [0u8; 348];

    // sizeof_hdr
    header[0..4].copy_from_slice(&348i32.to_le_bytes());
    // dim: [3, cols, rows, slices, 1, 1, 1, 1]
    header[40..42].copy_from_slice(&3i16.to_le_bytes());
    header[42..44].copy_from_slice(&cols.to_le_bytes());
    header[44..46].copy_from_slice(&rows.to_le_bytes());
    header[46..48].copy_from_slice(&slices.to_le_bytes());
    header[48..50].copy_from_slice(&1i16.to_le_bytes());
    header[50..52].copy_from_slice(&1i16.to_le_bytes());
    header[52..54].copy_from_slice(&1i16.to_le_bytes());
    header[54..56].copy_from_slice(&1i16.to_le_bytes());

    // datatype = 4 (INT16), bitpix = 16
    header[70..72].copy_from_slice(&4i16.to_le_bytes());
    header[72..74].copy_from_slice(&16i16.to_le_bytes());

    // pixdim
    header[76..80].copy_from_slice(&1.0f32.to_le_bytes());
    header[80..84].copy_from_slice(&pixdim_x.to_le_bytes());
    header[84..88].copy_from_slice(&pixdim_y.to_le_bytes());
    header[88..92].copy_from_slice(&pixdim_z.to_le_bytes());

    // vox_offset
    header[108..112].copy_from_slice(&352.0f32.to_le_bytes());

    // scl_slope = 1.0, scl_inter = 0.0
    header[112..116].copy_from_slice(&1.0f32.to_le_bytes());
    header[116..120].copy_from_slice(&0.0f32.to_le_bytes());

    // xyzt_units: mm + sec
    header[123] = 2 | 8;

    // sform_code = 1 (Scanner Anat)
    header[252..254].copy_from_slice(&1i16.to_le_bytes());
    // qform_code = 1
    header[254..256].copy_from_slice(&1i16.to_le_bytes());

    // srow_x, srow_y, srow_z (float32)
    let srow_offsets = [280, 296, 312];
    for (row_idx, &offset) in srow_offsets.iter().enumerate() {
        for col_idx in 0..4 {
            let val = affine[row_idx][col_idx] as f32;
            let pos = offset + col_idx * 4;
            header[pos..pos + 4].copy_from_slice(&val.to_le_bytes());
        }
    }

    // magic: "n+1\0"
    header[344..348].copy_from_slice(b"n+1\0");

    let file = std::fs::File::create(path)?;
    let mut gz = GzEncoder::new(file, Compression::fast());

    gz.write_all(&header)?;
    // 4-byte extension pad
    gz.write_all(&[0u8; 4])?;

    // Write voxel data as raw bytes
    let byte_slice = unsafe {
        std::slice::from_raw_parts(data.as_ptr() as *const u8, data.len() * 2)
    };
    gz.write_all(byte_slice)?;
    gz.finish()?;

    Ok(())
}
