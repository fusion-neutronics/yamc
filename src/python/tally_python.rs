#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use crate::tally::Tally;

#[cfg(feature = "pyo3")]
#[pyclass(name = "Tally")]
#[derive(Clone)]
pub struct PyTally {
    pub inner: Tally,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyTally {
    #[new]
    pub fn new() -> Self {
        PyTally {
            inner: Tally::new(),
        }
    }
    
    #[getter]
    pub fn score(&self) -> i32 {
        self.inner.score
    }
    
    #[setter]
    pub fn set_score(&mut self, score: i32) {
        self.inner.score = score;
    }
    
    #[getter]
    pub fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }
    
    #[setter]  
    pub fn set_name(&mut self, name: Option<String>) {
        self.inner.name = name;
    }
    
    #[getter]
    pub fn id(&self) -> Option<u32> {
        self.inner.id
    }
    
    #[setter]
    pub fn set_id(&mut self, id: Option<u32>) {
        self.inner.id = id;
    }

    #[getter]
    pub fn units(&self) -> String {
        self.inner.units.clone()
    }
    
    #[getter]
    pub fn mean(&self) -> f64 {
        self.inner.mean
    }
    
    #[getter]
    pub fn std_dev(&self) -> f64 {
        self.inner.std_dev
    }
    
    #[getter]
    pub fn rel_error(&self) -> f64 {
        self.inner.rel_error
    }
    
    #[getter]
    pub fn n_batches(&self) -> u32 {
        self.inner.n_batches
    }
    
    #[getter]
    pub fn particles_per_batch(&self) -> u32 {
        self.inner.particles_per_batch
    }
    
    #[getter]
    pub fn batch_data(&self) -> Vec<u32> {
        self.inner.batch_data.clone()
    }
    
    pub fn total_count(&self) -> u32 {
        self.inner.total_count()
    }
    
    fn __repr__(&self) -> String {
        format!("Tally(score={}, name={:?}, id={:?})", self.inner.score, self.inner.name, self.inner.id)
    }
    
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[cfg(feature = "pyo3")]
impl From<Tally> for PyTally {
    fn from(tally: Tally) -> Self {
        PyTally { inner: tally }
    }
}

#[cfg(feature = "pyo3")]
pub fn register_tally_classes(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTally>()?;
    Ok(())
}
