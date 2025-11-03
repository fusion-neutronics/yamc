use crate::stats::AngularDistribution;
use pyo3::prelude::*;

#[pyclass(name = "Isotropic")]
#[derive(Clone)]
pub struct PyIsotropic {
    pub inner: AngularDistribution,
}

#[pymethods]
impl PyIsotropic {
    #[new]
    pub fn new() -> Self {
        PyIsotropic {
            inner: AngularDistribution::Isotropic,
        }
    }

    pub fn sample(&self) -> [f64; 3] {
        // Python bindings use thread_rng for backwards compatibility
        let mut rng = rand::thread_rng();
        self.inner.sample(&mut rng)
    }

    pub fn __repr__(&self) -> String {
        "Isotropic()".to_string()
    }
}

#[pyclass(name = "Monodirectional")]
#[derive(Clone)]
pub struct PyMonodirectional {
    pub inner: AngularDistribution,
}

#[pymethods]
impl PyMonodirectional {
    #[new]
    pub fn new(reference_uvw: [f64; 3]) -> Self {
        PyMonodirectional {
            inner: AngularDistribution::new_monodirectional(
                reference_uvw[0],
                reference_uvw[1],
                reference_uvw[2],
            ),
        }
    }

    pub fn sample(&self) -> [f64; 3] {
        // Python bindings use thread_rng for backwards compatibility
        let mut rng = rand::thread_rng();
        self.inner.sample(&mut rng)
    }

    #[getter]
    pub fn reference_uvw(&self) -> [f64; 3] {
        match &self.inner {
            AngularDistribution::Monodirectional { reference_uvw } => *reference_uvw,
            AngularDistribution::Isotropic => panic!("Cannot get reference_uvw from Isotropic"),
        }
    }

    #[setter]
    pub fn set_reference_uvw(&mut self, reference_uvw: [f64; 3]) {
        self.inner = AngularDistribution::new_monodirectional(
            reference_uvw[0],
            reference_uvw[1],
            reference_uvw[2],
        );
    }

    pub fn __repr__(&self) -> String {
        match &self.inner {
            AngularDistribution::Monodirectional { reference_uvw } => {
                format!("Monodirectional(reference_uvw={:?})", reference_uvw)
            }
            AngularDistribution::Isotropic => {
                panic!("Invalid state: Monodirectional wrapper contains Isotropic")
            }
        }
    }
}

pub fn register_stats_module(py: Python, parent_module: &PyModule) -> PyResult<()> {
    let stats_module = PyModule::new(py, "stats")?;
    stats_module.add_class::<PyIsotropic>()?;
    stats_module.add_class::<PyMonodirectional>()?;
    parent_module.add_submodule(stats_module)?;
    Ok(())
}
