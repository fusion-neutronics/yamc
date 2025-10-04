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
    
    fn __repr__(&self) -> String {
        format!("Tally(score={}, name={:?}, id={:?})", self.inner.score, self.inner.name, self.inner.id)
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
