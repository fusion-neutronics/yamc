use crate::python::source_python::PyIndependentSource;
use crate::settings::Settings;
use pyo3::prelude::*;

#[pyclass(name = "Settings")]
#[derive(Clone)]
pub struct PySettings {
    pub inner: Settings,
}

#[pymethods]
impl PySettings {
    #[new]
    pub fn new(particles: usize, batches: usize, source: PyIndependentSource) -> Self {
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
    pub fn source(&self) -> PyIndependentSource {
        PyIndependentSource {
            inner: self.inner.source.clone(),
        }
    }

    #[setter]
    pub fn set_source(&mut self, value: PyIndependentSource) {
        self.inner.source = value.inner.clone();
    }
}
