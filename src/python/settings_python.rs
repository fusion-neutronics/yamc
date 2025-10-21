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
    #[pyo3(signature = (particles, batches, source, seed=None))]
    pub fn new(particles: usize, batches: usize, source: PyIndependentSource, seed: Option<u64>) -> Self {
        PySettings {
            inner: Settings {
                particles,
                batches,
                source: source.inner.clone(),
                seed,
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

    #[getter]
    pub fn seed(&self) -> Option<u64> {
        self.inner.seed
    }

    #[setter]
    pub fn set_seed(&mut self, value: Option<u64>) {
        self.inner.seed = value;
    }
}
