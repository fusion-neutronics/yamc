use crate::tallies::{CellFilter, EnergyFilter, MaterialFilter};
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
            Ok(filter) => Ok(PyMaterialFilter { internal: filter }),
            Err(_) => Err(pyo3::exceptions::PyValueError::new_err(
                "Cannot create MaterialFilter for material with no ID - assign a material_id first",
            )),
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

#[pyclass(name = "EnergyFilter")]
#[derive(Clone)]
pub struct PyEnergyFilter {
    pub internal: EnergyFilter,
}

#[pymethods]
impl PyEnergyFilter {
    /// Create a new EnergyFilter with energy bin boundaries
    ///
    /// Args:
    ///     bins (list[float]): Energy bin boundaries in eV, must be in ascending order
    ///
    /// Returns:
    ///     EnergyFilter: A new EnergyFilter with the specified energy bins
    ///
    /// Example:
    ///     # Create 3 energy bins: [0-1 MeV), [1-10 MeV), [10-20 MeV]
    ///     energy_filter = mc.EnergyFilter([0.0, 1e6, 10e6, 20e6])
    #[new]
    fn new(bins: Vec<f64>) -> PyResult<Self> {
        // Validate manually before calling EnergyFilter::new to avoid panic
        if bins.len() < 2 {
            return Err(pyo3::exceptions::PyValueError::new_err(
                "EnergyFilter requires at least 2 bin boundaries (to create at least 1 bin)",
            ));
        }
        
        // Verify bins are in ascending order
        for i in 1..bins.len() {
            if bins[i] <= bins[i - 1] {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "Energy bins must be in strictly ascending order",
                ));
            }
        }
        
        Ok(PyEnergyFilter {
            internal: EnergyFilter::new(bins),
        })
    }

    /// Get the energy bin boundaries
    #[getter]
    fn bins(&self) -> Vec<f64> {
        self.internal.bins.clone()
    }

    /// Get the number of energy bins (number of boundaries - 1)
    fn num_bins(&self) -> usize {
        self.internal.num_bins()
    }

    /// String representation of the EnergyFilter
    fn __repr__(&self) -> String {
        format!("EnergyFilter(bins={:?})", self.internal.bins)
    }

    /// String representation of the EnergyFilter
    fn __str__(&self) -> String {
        format!("EnergyFilter with {} bins", self.internal.num_bins())
    }
}

