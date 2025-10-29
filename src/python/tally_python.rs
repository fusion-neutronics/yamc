#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use pyo3::types::PyAny;
#[cfg(feature = "pyo3")]
use crate::tally::Tally;
#[cfg(feature = "pyo3")]
use std::sync::Arc;

#[cfg(feature = "pyo3")]
#[pyclass(name = "Tally", unsendable)]
#[derive(Clone)]
pub struct PyTally {
    pub inner: Arc<Tally>,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyTally {
    #[new]
    pub fn new() -> Self {
        PyTally {
            inner: Arc::new(Tally::new()),
        }
    }
    
    #[getter]
    pub fn scores(&self) -> Vec<i32> {
        self.inner.scores.clone()
    }

    #[setter]
    pub fn set_scores(&mut self, scores: Vec<i32>) {
        if let Some(tally) = Arc::get_mut(&mut self.inner) {
            tally.scores = scores;
        } else {
            panic!("Cannot modify tally: multiple references exist");
        }
    }

    #[getter]
    pub fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }

    #[setter]
    pub fn set_name(&mut self, name: Option<String>) {
        if let Some(tally) = Arc::get_mut(&mut self.inner) {
            tally.name = name;
        } else {
            panic!("Cannot modify tally: multiple references exist");
        }
    }

    #[getter]
    pub fn id(&self) -> Option<u32> {
        self.inner.id
    }

    #[setter]
    pub fn set_id(&mut self, id: Option<u32>) {
        if let Some(tally) = Arc::get_mut(&mut self.inner) {
            tally.id = id;
        } else {
            panic!("Cannot modify tally: multiple references exist");
        }
    }

    #[getter]
    pub fn units(&self) -> String {
        self.inner.units.clone()
    }
    
    #[getter]
    pub fn mean(&self) -> f64 {
        self.inner.get_mean()
    }

    #[getter]
    pub fn std_dev(&self) -> f64 {
        self.inner.get_std_dev()
    }

    #[getter]
    pub fn rel_error(&self) -> f64 {
        self.inner.get_rel_error()
    }

    #[getter]
    pub fn n_batches(&self) -> u32 {
        use std::sync::atomic::Ordering;
        self.inner.n_batches.load(Ordering::Relaxed)
    }

    #[getter]
    pub fn particles_per_batch(&self) -> u32 {
        use std::sync::atomic::Ordering;
        self.inner.particles_per_batch.load(Ordering::Relaxed)
    }

    #[getter]
    pub fn batch_data(&self) -> Vec<u64> {
        use std::sync::atomic::Ordering;
        self.inner.batch_data.lock().unwrap()
            .iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect()
    }
    
    #[getter]
    pub fn filters(&self) -> Vec<PyObject> {
        use pyo3::Python;
        use crate::python::filter_python::rust_filter_to_py;
        Python::with_gil(|py| {
            self.inner.filters
                .iter()
                .map(|f| rust_filter_to_py(py, f))
                .collect()
        })
    }
    
    #[setter]
    pub fn set_filters(&mut self, filters: Vec<&PyAny>) -> PyResult<()> {
        use crate::python::filter_python::py_filter_to_rust;
        let mut filter_objects = Vec::new();
        for filter_obj in filters {
            filter_objects.push(py_filter_to_rust(filter_obj)?);
        }
        if let Some(tally) = Arc::get_mut(&mut self.inner) {
            tally.filters = filter_objects;
            if let Err(err) = tally.validate() {
                return Err(pyo3::exceptions::PyValueError::new_err(err));
            }
        } else {
            return Err(pyo3::exceptions::PyRuntimeError::new_err(
                "Cannot modify tally: multiple references exist"
            ));
        }
        Ok(())
    }
    
    pub fn total_count(&self) -> u64 {
        self.inner.total_count()
    }
    
    fn __repr__(&self) -> String {
    format!("Tally(scores={:?}, name={:?}, id={:?})", self.inner.scores, self.inner.name, self.inner.id)
    }
    
    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[cfg(feature = "pyo3")]
impl From<Tally> for PyTally {
    fn from(tally: Tally) -> Self {
        PyTally { inner: Arc::new(tally) }
    }
}

#[cfg(feature = "pyo3")]
pub fn register_tally_classes(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTally>()?;
    Ok(())
}
