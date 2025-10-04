use crate::model::Model;
use crate::python::geometry_python::PyGeometry;
use crate::python::settings_python::PySettings;
use crate::python::source_python::PySource;
use crate::python::tally_python::PyTally;
use crate::tally::Tally;
use pyo3::prelude::*;

#[pyclass(name = "Model")]
#[derive(Clone)]
pub struct PyModel {
    pub inner: Model,
}

#[pymethods]
impl PyModel {
    #[new]
    #[pyo3(signature = (geometry, settings, tally=None))]
    pub fn new(geometry: PyGeometry, settings: PySettings, tally: Option<Vec<PyTally>>) -> Self {
        let tallies = if let Some(py_tallies) = tally {
            py_tallies.into_iter().map(|py_tally| py_tally.inner).collect()
        } else {
            Vec::new()
        };
        
        PyModel {
            inner: Model {
                geometry: geometry.inner.clone(),
                settings: settings.inner.clone(),
                tallies,
            },
        }
    }
    #[getter]
    pub fn geometry(&self) -> PyGeometry {
        PyGeometry {
            inner: self.inner.geometry.clone(),
        }
    }

    #[getter]
    pub fn settings(&self) -> PySettings {
        PySettings {
            inner: self.inner.settings.clone(),
        }
    }

    pub fn run(&self) -> Vec<PyTally> {
        let tallies = self.inner.run();
        tallies.into_iter().map(|t| t.into()).collect()
    }
}
