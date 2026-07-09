use dicom_core::Tag;
use dicom_core::header::Header;
use dicom_object::open_file;

fn main() {
    let path = std::env::args().nth(1).expect("Usage: debug_bval <dicom_file>");
    let obj = open_file(&path).expect("Failed to open DICOM file");
    
    // Tag (0018,9087) = DiffusionBValue
    let tag = Tag(0x0018, 0x9087);
    
    match obj.element(tag) {
        Ok(elem) => {
            println!("Tag found: {:?}", elem.tag());
            println!("VR: {:?}", elem.vr());
            println!("Value: {:?}", elem.value());
            
            // Try float64
            match elem.to_float64() {
                Ok(v) => println!("to_float64: {}", v),
                Err(e) => println!("to_float64 FAILED: {}", e),
            }
            
            // Try int
            match elem.to_int::<i32>() {
                Ok(v) => println!("to_int<i32>: {}", v),
                Err(e) => println!("to_int FAILED: {}", e),
            }
            
            // Try to_str
            match elem.to_str() {
                Ok(v) => println!("to_str: {:?}", v),
                Err(e) => println!("to_str FAILED: {}", e),
            }
            
            // Try multi_float64
            match elem.to_multi_float64() {
                Ok(v) => println!("to_multi_float64: {:?}", v),
                Err(e) => println!("to_multi_float64 FAILED: {}", e),
            }
        }
        Err(e) => {
            println!("Tag NOT found: {}", e);
            
            // Try alternative b-value tags
            for (g, e_tag, name) in &[
                (0x0043u16, 0x1039u16, "GE private"),
                (0x0019, 0x100C, "Siemens private"),
            ] {
                match obj.element(Tag(*g, *e_tag)) {
                    Ok(elem) => println!("  {} ({:04x},{:04x}): {:?}", name, g, e_tag, elem.value()),
                    Err(_) => println!("  {} ({:04x},{:04x}): NOT found", name, g, e_tag),
                }
            }
        }
    }
    
    // Also print manufacturer
    if let Ok(elem) = obj.element(Tag(0x0008, 0x0070)) {
        println!("Manufacturer: {:?}", elem.to_str().unwrap_or_default());
    }
    if let Ok(elem) = obj.element(Tag(0x0008, 0x103E)) {
        println!("SeriesDesc: {:?}", elem.to_str().unwrap_or_default());
    }
}
