use crate::source::Source;
use pyo3::prelude::*;

#[pyclass(name = "Source")]
#[derive(Clone)]
pub struct PySource {
    pub inner: Source,
}

#[pymethods]
impl PySource {
    #[new]
    pub fn new(position: [f64; 3], direction: [f64; 3], energy: f64) -> Self {
        PySource {
            inner: Source {
                position,
                direction,
                energy,
            },
        }
    }

    pub fn sample(&self) -> crate::python::particle_python::PyParticle {
        crate::python::particle_python::PyParticle {
            inner: self.inner.sample(),
        }
    }

    #[getter]
    pub fn position(&self) -> [f64; 3] {
        self.inner.position
    }
    #[getter]
    pub fn direction(&self) -> [f64; 3] {
        self.inner.direction
    }
    #[getter]
    pub fn energy(&self) -> f64 {
        self.inner.energy
    }
}
