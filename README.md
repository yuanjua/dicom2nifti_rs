# dicom2nifti-rs

A high-performance, parallel command-line tool written in Rust to convert DICOM image series into NIfTI-1 format (`.nii.gz`). It features robust slice grouping, orientation alignment, spatial transformation (LPS to RAS), subvolume/phase/b-value extraction, and parallel processing.

## Features

- **Parallel Processing**: Utilizes [Rayon](https://github.com/rayon-rs/rayon) to process multiple slice groups concurrently, leveraging all available CPU cores.
- **Robust Grouping & Sorting**: 
  - Groups DICOM slices by directory, Series Instance UID, and subvolume attributes.
  - Sorts slices along their calculated normal vector to guarantee correct 3D volume reconstruction.
- **Multi-Volume & Phase Splitting**:
  - **DWI (Diffusion Weighted Imaging)**: Auto-detects b-values using standard `DiffusionBValue` (0018,9087) as well as Siemens (0019,100c) and GE (0043,1039) private tags. Splitting is performed automatically into subvolumes (e.g. `b1000`, `b2000`).
  - **DCE / Dynamic Series**: Recognizes temporal acquisition phases via `TemporalPositionIdentifier` (0020,0100) and splits them (e.g., `phase1`, `phase2`).
- **Spatial Alignment (LPS to RAS)**:
  - Reconstructs 3D/4D spatial affines.
  - Corrects coordinate spaces by transforming DICOM LPS (Left, Posterior, Superior) to NIfTI RAS (Right, Anterior, Superior).
- **Inconsistency Diagnostics & Splitting**:
  - Detects orientation variations within a series and splits them into distinct plane groups (e.g. `plane0`, `plane1`).
  - Detects non-uniform slice spacing and splits groups if gaps exceed five times the median gap (e.g. `part0`, `part1`).
- **Detailed Audit Logs**:
  - **`conversion_mapping.csv`**: Contains complete metadata (Patient ID, Study Description, Series Number, Manufacturer, Spatial Info, Spacing, Shape, and File Count) of all successfully converted volumes.
  - **`error_log.csv`**: Logs any files or groups that failed to convert along with their errors.

---

## Installation

Ensure you have Rust and Cargo installed (see [rustup.rs](https://rustup.rs/)).

1. Clone or download the repository.
2. Build the project in release mode:
   ```bash
   cargo build --release
   ```
   The compiled binary will be available at [target/release/dicom2nifti-rs](file:///home/caiyuanzhou/tmp/dicom2nifti_rs/target/release/dicom2nifti-rs).

---

## Usage

```bash
cargo run --release -- <INPUT_DIR> -o <OUTPUT_DIR> [FLAGS]
```

### CLI Arguments & Options

| Argument/Option | Short | Description |
| :--- | :--- | :--- |
| `<INPUT_DIR>` | | Input directory containing DICOM files (searched recursively). |
| `-o, --output` | `-o` | Output directory where `.nii.gz` and CSV logs will be saved. |
| `-j, --threads` | `-j` | Number of parallel threads to use. Defaults to `0` (auto-detects based on system cores). |
| `--help` | `-h` | Prints help information. |
| `--version` | `-V` | Prints version information. |

### Example Execution

```bash
cargo run --release -- /path/to/dicom/dataset -o /path/to/output_nifti
```

---

## Codebase Tour

- **[src/main.rs](file:///home/caiyuanzhou/tmp/dicom2nifti_rs/src/main.rs)**: The CLI entry point, WalkDir-based file search, header reading, grouping logic, and Rayon parallel orchestration.
- **[src/spatial.rs](file:///home/caiyuanzhou/tmp/dicom2nifti_rs/src/spatial.rs)**: Contains coordinate transformation logic and computes the 4x4 affine mapping matrix for LPS → RAS.
- **[src/nifti_write.rs](file:///home/caiyuanzhou/tmp/dicom2nifti_rs/src/nifti_write.rs)**: Implements writing raw volume bytes to a `.nii.gz` file with a compliant NIfTI-1 header.
- **[examples/](file:///home/caiyuanzhou/tmp/dicom2nifti_rs/examples)**:
  - **[debug_bval.rs](file:///home/caiyuanzhou/tmp/dicom2nifti_rs/examples/debug_bval.rs)**: A helper tool to open a single DICOM file and print its parsed b-value tags (useful for diagnosing vendor-specific tags).
  - **[debug_group.rs](file:///home/caiyuanzhou/tmp/dicom2nifti_rs/examples/debug_group.rs)**: A diagnostic tool to scan a directory and output the b-value distribution of DICOM slices.

To run the examples:
```bash
cargo run --example debug_bval -- /path/to/file.dcm
cargo run --example debug_group -- /path/to/series_directory
```
