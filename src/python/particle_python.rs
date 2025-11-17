use crate::particle::Particle;
use pyo3::prelude::*;

#[pyclass(name = "Particle")]
#[derive(Clone)]
pub struct PyParticle {
    pub inner: Particle,
}

#[pymethods]
impl PyParticle {
    #[new]
    pub fn new(position: [f64; 3], direction: [f64; 3], energy: f64, alive: Option<bool>, id: Option<usize>) -> Self {
        PyParticle {
            inner: Particle {
                position,
                direction,
                energy,
                alive: alive.unwrap_or(true),
                id: id.unwrap_or(0),
                current_cell_index: None, // Will be set during transport
            },
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
    #[getter]
    pub fn alive(&self) -> bool {
        self.inner.alive
    }
    #[getter]
    pub fn id(&self) -> usize {
        self.inner.id
    }
}
