use pyo3::prelude::*;
use crate::model::Model;
use crate::python::geometry_python::PyGeometry;
use crate::python::source_python::PySource;
use crate::python::settings_python::PySettings;

#[pyclass(name = "Model")]
#[derive(Clone)]
pub struct PyModel {
    pub inner: Model,
}

#[pymethods]
impl PyModel {
    #[new]
    pub fn new(geometry: PyGeometry, settings: PySettings) -> Self {
        PyModel {
            inner: Model {
                geometry: geometry.inner.clone(),
                settings: settings.inner.clone(),
            },
        }
    }
    #[getter]
    pub fn geometry(&self) -> PyGeometry {
        PyGeometry { inner: self.inner.geometry.clone() }
    }

    #[getter]
    pub fn settings(&self) -> PySettings {
        PySettings { inner: self.inner.settings.clone() }
    }

    pub fn run(&self) {
        self.inner.run();
    }
}
