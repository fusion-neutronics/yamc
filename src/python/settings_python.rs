use pyo3::prelude::*;
use crate::settings::Settings;
use crate::python::source_python::PySource;

#[pyclass(name = "Settings")]
#[derive(Clone)]
pub struct PySettings {
    pub inner: Settings,
}

#[pymethods]
impl PySettings {
    #[new]
    pub fn new(particles: usize, batches: usize, source: PySource) -> Self {
        PySettings {
            inner: Settings {
                particles,
                batches,
                source: source.inner.clone(),
            },
        }
    }
    #[getter]
    pub fn particles(&self) -> usize {
        self.inner.particles
    }
    #[getter]
    pub fn batches(&self) -> usize {
        self.inner.batches
    }
    #[getter]
    pub fn source(&self) -> PySource {
        PySource { inner: self.inner.source.clone() }
    }

    #[setter]
    pub fn set_source(&mut self, value: PySource) {
        self.inner.source = value.inner.clone();
    }
}
