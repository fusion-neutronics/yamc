#[cfg(feature = "pyo3")]
use crate::tally::{Tally, Score};
#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use pyo3::types::PyAny;
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
    pub fn scores(&self) -> Vec<PyObject> {
        Python::with_gil(|py| {
            self.inner.scores.iter().map(|score| match score {
                Score::MT(mt) => mt.into_py(py),
                Score::Flux => "flux".into_py(py),
            }).collect()
        })
    }

    #[setter]
    pub fn set_scores(&mut self, scores: &PyAny) -> PyResult<()> {
        use std::sync::atomic::AtomicU64;
        use std::sync::{Arc as StdArc, Mutex};
        
        // Convert Python input to Vec<Score>
        let score_list: Vec<Score> = if let Ok(ints) = scores.extract::<Vec<i32>>() {
            ints.into_iter().map(Score::MT).collect()
        } else if let Ok(items) = scores.extract::<Vec<&PyAny>>() {
            items.iter().map(|item| {
                if let Ok(i) = item.extract::<i32>() {
                    Ok(Score::MT(i))
                } else if let Ok(s) = item.extract::<String>() {
                    if s == "flux" {
                        Ok(Score::Flux)
                    } else {
                        Err(PyErr::new::<pyo3::exceptions::PyValueError, _>(format!("Unknown score: {}", s)))
                    }
                } else {
                    Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>("Score must be int or string"))
                }
            }).collect::<PyResult<Vec<Score>>>()?
        } else {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>("scores must be a list"));
        };

        let tally = Arc::get_mut(&mut self.inner)
            .ok_or_else(|| PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(
                "Cannot modify tally: multiple references exist"
            ))?;
        
        // Set scores and pre-allocate storage accounting for energy bins
        tally.scores = score_list;
        let num_bins = tally.scores.len() * tally.num_energy_bins();
        
        tally.batch_data = (0..num_bins)
            .map(|_| Mutex::new(StdArc::new(Vec::new())))
            .collect();
        tally.mean = (0..num_bins)
            .map(|_| AtomicU64::new(0))
            .collect();
        tally.std_dev = (0..num_bins)
            .map(|_| AtomicU64::new(0))
            .collect();
        tally.rel_error = (0..num_bins)
            .map(|_| AtomicU64::new(0))
            .collect();
        
        Ok(())
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
    pub fn mean(&self) -> Vec<f64> {
        self.inner.get_mean()
    }

    #[getter]
    pub fn std_dev(&self) -> Vec<f64> {
        self.inner.get_std_dev()
    }

    #[getter]
    pub fn rel_error(&self) -> Vec<f64> {
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
    pub fn batch_data(&self) -> Vec<Vec<u64>> {
        use std::sync::atomic::Ordering;
        self.inner
            .batch_data
            .iter()
            .map(|batch_mutex| {
                batch_mutex
                    .lock()
                    .unwrap()
                    .iter()
                    .map(|a| a.load(Ordering::Relaxed))
                    .collect()
            })
            .collect()
    }

    #[getter]
    pub fn filters(&self) -> Vec<PyObject> {
        use crate::python::filter_python::rust_filter_to_py;
        use pyo3::Python;
        Python::with_gil(|py| {
            self.inner
                .filters
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
                "Cannot modify tally: multiple references exist",
            ));
        }
        Ok(())
    }

    pub fn total_count(&self) -> Vec<u64> {
        self.inner.total_count()
    }

    fn __repr__(&self) -> String {
        // Format scores as Python list: integers for MT, "flux" for Flux
        let scores_str = self.inner.scores.iter().map(|score| match score {
            Score::MT(mt) => mt.to_string(),
            Score::Flux => "\"flux\"".to_string(),
        }).collect::<Vec<_>>().join(", ");
        
        format!(
            "Tally(scores=[{}], name={:?}, id={:?})",
            scores_str, self.inner.name, self.inner.id
        )
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}

#[cfg(feature = "pyo3")]
impl From<Tally> for PyTally {
    fn from(tally: Tally) -> Self {
        PyTally {
            inner: Arc::new(tally),
        }
    }
}

#[cfg(feature = "pyo3")]
pub fn register_tally_classes(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<PyTally>()?;
    Ok(())
}
