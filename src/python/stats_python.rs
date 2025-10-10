use crate::stats::{AngularDistribution, Isotropic, Monodirectional};
use pyo3::prelude::*;

#[pyclass(name = "Isotropic")]
#[derive(Clone)]
pub struct PyIsotropic {
    pub inner: Isotropic,
}

#[pymethods]
impl PyIsotropic {
    #[new]
    pub fn new() -> Self {
        PyIsotropic {
            inner: Isotropic,
        }
    }

    pub fn sample(&self) -> [f64; 3] {
        self.inner.sample()
    }

    pub fn __repr__(&self) -> String {
        "Isotropic()".to_string()
    }
}

#[pyclass(name = "Monodirectional")]
#[derive(Clone)]
pub struct PyMonodirectional {
    pub inner: Monodirectional,
}

#[pymethods]
impl PyMonodirectional {
    #[new]
    pub fn new(reference_uvw: [f64; 3]) -> Self {
        PyMonodirectional {
            inner: Monodirectional::new(reference_uvw[0], reference_uvw[1], reference_uvw[2]),
        }
    }

    pub fn sample(&self) -> [f64; 3] {
        self.inner.sample()
    }

    #[getter]
    pub fn reference_uvw(&self) -> [f64; 3] {
        self.inner.reference_uvw
    }

    #[setter]
    pub fn set_reference_uvw(&mut self, reference_uvw: [f64; 3]) {
        self.inner.reference_uvw = reference_uvw;
    }

    pub fn __repr__(&self) -> String {
        format!("Monodirectional(reference_uvw={:?})", self.inner.reference_uvw)
    }
}

pub fn register_stats_module(py: Python, parent_module: &PyModule) -> PyResult<()> {
    let stats_module = PyModule::new(py, "stats")?;
    stats_module.add_class::<PyIsotropic>()?;
    stats_module.add_class::<PyMonodirectional>()?;
    parent_module.add_submodule(stats_module)?;
    Ok(())
}