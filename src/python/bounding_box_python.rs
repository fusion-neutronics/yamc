use pyo3::prelude::*;

#[pyclass(name = "BoundingBox")]
#[derive(Clone)]
pub struct PyBoundingBox {
    pub lower_left: [f64; 3],
    pub upper_right: [f64; 3],
}

#[pymethods]
impl PyBoundingBox {
    #[new]
    pub fn new(lower_left: [f64; 3], upper_right: [f64; 3]) -> Self {
        PyBoundingBox { lower_left, upper_right }
    }

    #[getter]
    pub fn lower_left(&self) -> [f64; 3] {
        self.lower_left
    }

    #[getter]
    pub fn upper_right(&self) -> [f64; 3] {
        self.upper_right
    }

    pub fn __repr__(&self) -> String {
        format!("BoundingBox(lower_left={:?}, upper_right={:?})", self.lower_left, self.upper_right)
    }
}
