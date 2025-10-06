use crate::filters::{CellFilter, MaterialFilter};
use crate::python::cell_python::PyCell;
use crate::python::material_python::PyMaterial;
use pyo3::prelude::*;

#[pyclass(name = "CellFilter")]
#[derive(Clone)]
pub struct PyCellFilter {
    pub internal: CellFilter,
}

#[pymethods]
impl PyCellFilter {
    /// Create a new CellFilter from a Cell object
    ///
    /// Args:
    ///     cell (Cell): The cell object to create the filter from
    ///
    /// Returns:
    ///     CellFilter: A new CellFilter that will match events in this cell
    #[new]
    fn new(cell: &PyCell) -> Self {
        // Extract the internal Cell from PyCell
        PyCellFilter {
            internal: CellFilter::new(&cell.inner),
        }
    }



    /// String representation of the CellFilter
    fn __repr__(&self) -> String {
        format!("CellFilter(cell_id={})", self.internal.cell_id)
    }

    /// String representation of the CellFilter
    fn __str__(&self) -> String {
        format!("CellFilter for cell {}", self.internal.cell_id)
    }
}

#[pyclass(name = "MaterialFilter")]
#[derive(Clone)]
pub struct PyMaterialFilter {
    pub internal: MaterialFilter,
}

#[pymethods]
impl PyMaterialFilter {
    /// Create a new MaterialFilter from a Material object
    ///
    /// Args:
    ///     material (Material): The material object to create the filter from
    ///
    /// Returns:
    ///     MaterialFilter: A new MaterialFilter that will match events in this material
    #[new]
    fn new(material: &PyMaterial) -> PyResult<Self> {
        match std::panic::catch_unwind(|| MaterialFilter::new(&material.internal)) {
            Ok(filter) => Ok(PyMaterialFilter {
                internal: filter,
            }),
            Err(_) => Err(pyo3::exceptions::PyValueError::new_err(
                "Cannot create MaterialFilter for material with no ID - assign a material_id first"
            ))
        }
    }

    /// String representation of the MaterialFilter
    fn __repr__(&self) -> String {
        format!("MaterialFilter(material_id={})", self.internal.material_id)
    }

    /// String representation of the MaterialFilter
    fn __str__(&self) -> String {
        format!("MaterialFilter for material {}", self.internal.material_id)
    }
}