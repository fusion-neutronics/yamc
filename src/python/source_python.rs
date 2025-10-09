use crate::source::IndependentSource;
use pyo3::prelude::*;

#[pyclass(name = "IndependentSource")]
#[derive(Clone)]
pub struct PyIndependentSource {
    pub inner: IndependentSource,
}

#[pymethods]
impl PyIndependentSource {
    #[new]
    pub fn new(space: [f64; 3], angle: [f64; 3], energy: f64) -> Self {
        PyIndependentSource {
            inner: IndependentSource {
                space,
                angle,
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
    pub fn space(&self) -> [f64; 3] {
        self.inner.space
    }
    #[getter]
    pub fn angle(&self) -> [f64; 3] {
        self.inner.angle
    }
    #[getter]
    pub fn energy(&self) -> f64 {
        self.inner.energy
    }
}
