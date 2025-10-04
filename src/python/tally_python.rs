#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use pyo3::exceptions::{PyIndexError, PyTypeError};
#[cfg(feature = "pyo3")]
use pyo3::types::PyAny;
#[cfg(feature = "pyo3")]
use crate::tally::{CountTally, Tally};

#[cfg(feature = "pyo3")]
#[pyclass]
#[derive(Clone)]
pub struct PyCountTally {
    #[pyo3(get)]
    pub name: String,
    #[pyo3(get)]
    pub units: String,
    #[pyo3(get)]
    pub batch_data: Vec<u32>,
    #[pyo3(get)]
    pub mean: f64,
    #[pyo3(get)]
    pub std_dev: f64,
    #[pyo3(get)]
    pub rel_error: f64,
    #[pyo3(get)]
    pub n_batches: u32,
    #[pyo3(get)]
    pub particles_per_batch: u32,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyCountTally {
    #[getter]
    pub fn total_count(&self) -> u32 {
        self.batch_data.iter().sum()
    }
    
    fn __str__(&self) -> String {
        format!(
            "Tally name: {}:\n  Mean: {:.6} per particle\n    Std Dev: {:.6} per particle\n    Rel Error: {:.4} ({:.2}%)\n    Batches: {}\n    Particles per batch: {}\n  Total {}: {}\n  Batch data: {:?}",
            self.name,
            self.mean,
            self.std_dev,
            self.rel_error,
            self.rel_error * 100.0,
            self.n_batches,
            self.particles_per_batch,
            self.units.to_lowercase(),
            self.total_count(),
            self.batch_data
        )
    }
    
    fn __repr__(&self) -> String {
        self.__str__()
    }
}

#[cfg(feature = "pyo3")]
impl From<CountTally> for PyCountTally {
    fn from(tally: CountTally) -> Self {
        PyCountTally {
            name: tally.name,
            units: tally.units,
            batch_data: tally.batch_data,
            mean: tally.mean,
            std_dev: tally.std_dev,
            rel_error: tally.rel_error,
            n_batches: tally.n_batches,
            particles_per_batch: tally.particles_per_batch,
        }
    }
}





#[cfg(feature = "pyo3")]
#[pyclass]
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
    pub fn score(&self) -> Vec<i32> {
        self.inner.score.clone()
    }
    
    #[setter]
    pub fn set_score(&mut self, score: &PyAny) -> PyResult<()> {
        // Try to convert to a single integer first
        if let Ok(single_score) = score.extract::<i32>() {
            self.inner.score = vec![single_score];
        }
        // If that fails, try to convert to a list of integers
        else if let Ok(score_list) = score.extract::<Vec<i32>>() {
            self.inner.score = score_list;
        }
        // If both fail, return an error
        else {
            return Err(PyTypeError::new_err("Score must be an integer or a list of integers"));
        }
        Ok(())
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
    
    fn __repr__(&self) -> String {
        format!("Tally(score={:?}, name={:?}, id={:?})", 
                self.inner.score, self.inner.name, self.inner.id)
    }
}

#[cfg(feature = "pyo3")]
pub fn register_tally_classes(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyCountTally>()?;
    // Register PyTally as "Tally" for user convenience
    m.add("Tally", _py.get_type::<PyTally>())?;
    Ok(())
}