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
    pub fn score(&self) -> i32 {
        self.inner.score
    }
    
    #[setter]
    pub fn set_score(&mut self, score: i32) {
        if let Some(tally) = Arc::get_mut(&mut self.inner) {
            tally.score = score;
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
        self.inner.mean.get()
    }

    #[getter]
    pub fn std_dev(&self) -> f64 {
        self.inner.std_dev.get()
    }

    #[getter]
    pub fn rel_error(&self) -> f64 {
        self.inner.rel_error.get()
    }

    #[getter]
    pub fn n_batches(&self) -> u32 {
        self.inner.n_batches.get()
    }
    
    #[getter]
    pub fn particles_per_batch(&self) -> u32 {
        self.inner.particles_per_batch.get()
    }
    
    #[getter]
    pub fn batch_data(&self) -> Vec<u64> {
        use std::sync::atomic::Ordering;
        self.inner.batch_data.borrow()
            .iter()
            .map(|a| a.load(Ordering::Relaxed))
            .collect()
    }
    
    #[getter]
    pub fn filters(&self) -> Vec<PyObject> {
        // Convert internal Filters to appropriate Python filter objects
        use pyo3::Python;
        Python::with_gil(|py| {
            self.inner.filters
                .iter()
                .map(|f| match f {
                    crate::tally::Filter::Cell(cell_filter) => {
                        let py_cell_filter = crate::python::filters_python::PyCellFilter {
                            internal: cell_filter.clone(),
                        };
                        py_cell_filter.into_py(py)
                    }
                    crate::tally::Filter::Material(material_filter) => {
                        let py_material_filter = crate::python::filters_python::PyMaterialFilter {
                            internal: material_filter.clone(),
                        };
                        py_material_filter.into_py(py)
                    }
                })
                .collect()
        })
    }
    
    #[setter]
    pub fn set_filters(&mut self, filters: Vec<&PyAny>) -> PyResult<()> {
        let mut filter_objects = Vec::new();
        
        for filter_obj in filters {
            if let Ok(cell_filter) = filter_obj.extract::<crate::python::filters_python::PyCellFilter>() {
                filter_objects.push(crate::tally::Filter::Cell(cell_filter.internal));
            } else if let Ok(material_filter) = filter_obj.extract::<crate::python::filters_python::PyMaterialFilter>() {
                filter_objects.push(crate::tally::Filter::Material(material_filter.internal));
            } else {
                return Err(pyo3::exceptions::PyTypeError::new_err(
                    "Filter must be either CellFilter or MaterialFilter"
                ));
            }
        }
        
        if let Some(tally) = Arc::get_mut(&mut self.inner) {
            tally.filters = filter_objects;

            // Validate the tally configuration
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
        format!("Tally(score={}, name={:?}, id={:?})", self.inner.score, self.inner.name, self.inner.id)
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
