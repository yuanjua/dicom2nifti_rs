use dicom_core::Tag;
use dicom_object::open_file;
use std::collections::HashMap;
use std::path::Path;

fn main() {
    let dir = std::env::args().nth(1).expect("Usage: debug_group <series_dir>");
    let dir = Path::new(&dir);
    
    let mut bval_counts: HashMap<i64, usize> = HashMap::new();
    
    for entry in std::fs::read_dir(dir).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        if !path.is_file() { continue; }
        let ext = path.extension().and_then(|s| s.to_str()).unwrap_or("");
        if !ext.eq_ignore_ascii_case("dcm") { continue; }
        
        let obj = match open_file(&path) {
            Ok(o) => o,
            Err(e) => {
                eprintln!("Failed to open {:?}: {}", path, e);
                continue;
            }
        };
        
        // Try reading b-value exactly as determine_subvolume does
        let mut bval: i64 = 0;
        let mut label = String::new();
        
        // Strategy 1: TemporalPositionIdentifier
        if let Ok(elem) = obj.element(Tag(0x0020, 0x0100)) {
            if let Ok(tp) = elem.to_int::<i32>() {
                if tp > 0 {
                    bval = tp as i64;
                    label = format!("phase{tp}");
                    println!("  {:?}: TP={} -> ({}, {})", path.file_name().unwrap(), tp, bval, label);
                    *bval_counts.entry(bval).or_insert(0) += 1;
                    continue;
                }
            }
        }
        
        // Strategy 2: DiffusionBValue (0018,9087)
        if let Ok(elem) = obj.element(Tag(0x0018, 0x9087)) {
            if let Ok(bv) = elem.to_float64() {
                bval = bv.round() as i64;
                label = format!("b{bval}");
                println!("  {:?}: FD bval={} -> ({}, {})", path.file_name().unwrap(), bv, bval, label);
                *bval_counts.entry(bval).or_insert(0) += 1;
                continue;
            }
            if let Ok(bv) = elem.to_int::<i32>() {
                bval = bv as i64;
                label = format!("b{bval}");
                println!("  {:?}: IS bval={} -> ({}, {})", path.file_name().unwrap(), bv, bval, label);
                *bval_counts.entry(bval).or_insert(0) += 1;
                continue;
            }
            eprintln!("  {:?}: tag exists but NEITHER float64 nor int worked! value={:?}", 
                      path.file_name().unwrap(), elem.value());
        }
        
        println!("  {:?}: no bval found -> (0, \"\")", path.file_name().unwrap());
        *bval_counts.entry(0).or_insert(0) += 1;
    }
    
    println!("\nB-value distribution:");
    for (bv, cnt) in &bval_counts {
        println!("  b={}: {} files", bv, cnt);
    }
}
