#[derive(Debug, Clone, PartialEq)]
pub struct BoundingBox {
    pub lower_left: [f64; 3],
    pub upper_right: [f64; 3],
    pub center: [f64; 3],
    pub width: [f64; 3],
}

impl BoundingBox {
    pub fn new(lower_left: [f64; 3], upper_right: [f64; 3]) -> Self {
        let center = [
            0.5 * (lower_left[0] + upper_right[0]),
            0.5 * (lower_left[1] + upper_right[1]),
            0.5 * (lower_left[2] + upper_right[2]),
        ];
        let width = [
            upper_right[0] - lower_left[0],
            upper_right[1] - lower_left[1],
            upper_right[2] - lower_left[2],
        ];
        BoundingBox {
            lower_left,
            upper_right,
            center,
            width,
        }
    }
}
