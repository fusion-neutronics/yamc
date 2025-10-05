use crate::filters::CellFilter;
use crate::cell::Cell;
use pyo3::prelude::*;
use pyo3::types::PyAny;

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

// Import the PyCell from the cell_python module
use crate::python::cell_python::PyCell;