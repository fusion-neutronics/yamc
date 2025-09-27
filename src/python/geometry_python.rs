use crate::geometry;
use crate::cell_python::PyCell;

#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
use pyo3::types::PySequence;

#[cfg_attr(feature = "pyo3", pyclass(name = "Geometry"))]
#[derive(Clone)]
pub struct PyGeometry {
    pub inner: geometry::Geometry,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyGeometry {
    #[new]
    pub fn new(cells: &pyo3::types::PyAny) -> pyo3::PyResult<Self> {
        let seq = cells.downcast::<PySequence>()?;
        let mut rust_cells = Vec::with_capacity(seq.len()?);
        for item in seq.iter()? {
            let pycell: PyCell = item?.extract()?;
            rust_cells.push(pycell.inner);
        }
        Ok(PyGeometry { inner: geometry::Geometry { cells: rust_cells } })
    }

    pub fn find_cell(&self, x: f64, y: f64, z: f64) -> Option<PyCell> {
        self.inner.find_cell((x, y, z)).cloned().map(|cell| PyCell { inner: cell })
    }
}
