use crate::model::Model;
use crate::python::geometry_python::PyGeometry;
use crate::python::settings_python::PySettings;
use crate::python::tally_python::PyTally;
use crate::tally::Tally;
use pyo3::prelude::*;

#[pyclass(name = "Model", unsendable)]
#[derive(Clone)]
pub struct PyModel {
    pub inner: Model,
}

#[pymethods]
impl PyModel {
    #[new]
    #[pyo3(signature = (geometry, settings, tallies=None))]
    pub fn new(geometry: PyGeometry, settings: PySettings, tallies: Option<Vec<PyTally>>) -> Self {
        let tallies = if let Some(py_tallies) = tallies {
            // Clone the Arc from each PyTally
            py_tallies.into_iter().map(|py_tally| py_tally.inner.clone()).collect()
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
    
    #[getter]
    pub fn tallies(&self) -> Vec<PyTally> {
        self.inner.tallies.iter().map(|t| PyTally { inner: t.clone() }).collect()
    }

    pub fn run(&mut self) {
        self.inner.run();
        // Tallies are updated in place - no return value
    }
}
