#[cfg(feature = "pyo3")]
use crate::filter::Filter;
#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use pyo3::types::PyAny;

#[cfg(feature = "pyo3")]
pub fn py_filter_to_rust(filter_obj: &PyAny) -> PyResult<Filter> {
    if let Ok(cell_filter) = filter_obj.extract::<crate::python::filters_python::PyCellFilter>() {
        Ok(Filter::Cell(cell_filter.internal))
    } else if let Ok(material_filter) =
        filter_obj.extract::<crate::python::filters_python::PyMaterialFilter>()
    {
        Ok(Filter::Material(material_filter.internal))
    } else if let Ok(energy_filter) =
        filter_obj.extract::<crate::python::filters_python::PyEnergyFilter>()
    {
        Ok(Filter::Energy(energy_filter.internal))
    } else {
        Err(pyo3::exceptions::PyTypeError::new_err(
            "Filter must be CellFilter, MaterialFilter, or EnergyFilter",
        ))
    }
}

#[cfg(feature = "pyo3")]
pub fn rust_filter_to_py(py: Python, filter: &Filter) -> PyObject {
    match filter {
        Filter::Cell(cell_filter) => {
            let py_cell_filter = crate::python::filters_python::PyCellFilter {
                internal: cell_filter.clone(),
            };
            py_cell_filter.into_py(py)
        }
        Filter::Material(material_filter) => {
            let py_material_filter = crate::python::filters_python::PyMaterialFilter {
                internal: material_filter.clone(),
            };
            py_material_filter.into_py(py)
        }
        Filter::Energy(energy_filter) => {
            let py_energy_filter = crate::python::filters_python::PyEnergyFilter {
                internal: energy_filter.clone(),
            };
            py_energy_filter.into_py(py)
        }
    }
}
