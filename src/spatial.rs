#[derive(Debug, Clone)]
pub struct SliceInfo {
    pub position: [f64; 3],
    pub orientation: [f64; 6],
    pub pixel_spacing: [f64; 2],
    pub slice_thickness: f64,
}

impl SliceInfo {
    pub fn slice_normal(&self) -> [f64; 3] {
        let o = &self.orientation;
        [
            o[1] * o[5] - o[2] * o[4],
            o[2] * o[3] - o[0] * o[5],
            o[0] * o[4] - o[1] * o[3],
        ]
    }

    pub fn dot_normal(&self, normal: &[f64; 3]) -> f64 {
        self.position[0] * normal[0]
            + self.position[1] * normal[1]
            + self.position[2] * normal[2]
    }
}

/// Compute the 4×4 affine matrix (DICOM LPS → NIfTI RAS) from first/last slice info.
pub fn compute_affine(
    first: &SliceInfo,
    last: &SliceInfo,
    rows: usize,
    cols: usize,
    n_slices: usize,
) -> [[f64; 4]; 4] {
    let o = &first.orientation;
    let ps = &first.pixel_spacing;

    // Row and column direction cosines (DICOM LPS)
    let row_x = o[0];
    let row_y = o[1];
    let row_z = o[2];
    let col_x = o[3];
    let col_y = o[4];
    let col_z = o[5];

    // Slice direction from first→last position
    let (slice_x, slice_y, slice_z, slice_spacing) = if n_slices > 1 {
        let dx = last.position[0] - first.position[0];
        let dy = last.position[1] - first.position[1];
        let dz = last.position[2] - first.position[2];
        let dist = (dx * dx + dy * dy + dz * dz).sqrt();
        let sp = dist / (n_slices - 1) as f64;
        if dist > 1e-10 {
            (dx / dist * sp, dy / dist * sp, dz / dist * sp, sp)
        } else {
            let n = first.slice_normal();
            let st = first.slice_thickness;
            (n[0] * st, n[1] * st, n[2] * st, st)
        }
    } else {
        let n = first.slice_normal();
        let st = first.slice_thickness;
        (n[0] * st, n[1] * st, n[2] * st, st)
    };

    // LPS→RAS: negate x and y
    let mut m = [[0.0f64; 4]; 4];
    m[0][0] = -row_x * ps[1];
    m[0][1] = -col_x * ps[0];
    m[0][2] = -slice_x;
    m[0][3] = -first.position[0];

    m[1][0] = -row_y * ps[1];
    m[1][1] = -col_y * ps[0];
    m[1][2] = -slice_y;
    m[1][3] = -first.position[1];

    m[2][0] = row_z * ps[1];
    m[2][1] = col_z * ps[0];
    m[2][2] = slice_z;
    m[2][3] = first.position[2];

    m[3][3] = 1.0;
    let _ = slice_spacing;
    m
}
