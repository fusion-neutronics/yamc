use crate::source::IndependentSource;
use pyo3::prelude::*;
use pyo3::types::PyAny;
use crate::python::stats_python::{PyIsotropic, PyMonodirectional};
use crate::stats::{AngularDistribution, Isotropic, Monodirectional};

impl<'source> FromPyObject<'source> for AngularDistribution {
    fn extract(ob: &'source PyAny) -> PyResult<Self> {
        if let Ok(mono) = ob.extract::<PyMonodirectional>() {
            Ok(AngularDistribution::Monodirectional { 
                reference_uvw: mono.inner.reference_uvw 
            })
        } else if let Ok(_iso) = ob.extract::<PyIsotropic>() {
            Ok(AngularDistribution::Isotropic)
        } else {
            Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "Expected PyIsotropic or PyMonodirectional"
            ))
        }
    }
}

impl IntoPy<PyObject> for AngularDistribution {
    fn into_py(self, py: Python<'_>) -> PyObject {
        match self {
            AngularDistribution::Isotropic => {
                PyIsotropic { inner: Isotropic }.into_py(py)
            },
            AngularDistribution::Monodirectional { reference_uvw } => {
                let mono = Monodirectional { reference_uvw };
                PyMonodirectional { inner: mono }.into_py(py)
            }
        }
    }
}

#[pyclass(name = "IndependentSource")]
#[derive(Clone)]
pub struct PyIndependentSource {
    pub inner: IndependentSource,
}

#[pymethods]
impl PyIndependentSource {
    #[new]
    #[pyo3(signature = (*, space=None, energy=None, angle=None))]
    pub fn new(
        space: Option<[f64; 3]>,
        energy: Option<f64>, 
        angle: Option<AngularDistribution>
    ) -> PyResult<Self> {
        let mut instance = Self { inner: IndependentSource::new() };
        
        if let Some(val) = space { instance.set_space(val); }
        if let Some(val) = energy { instance.set_energy(val); }
        if let Some(val) = angle { instance.set_angle_direct(val); }
        
        Ok(instance)
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

    #[setter(space)]
    pub fn set_space(&mut self, space: [f64; 3]) {
        self.inner.space = space;
    }

    #[getter]
    pub fn energy(&self) -> f64 {
        self.inner.energy
    }

    #[setter(energy)]
    pub fn set_energy(&mut self, energy: f64) {
        self.inner.energy = energy;
    }

    #[getter]
    pub fn angle(&self) -> AngularDistribution {
        self.inner.angle.clone()
    }

    fn set_angle_direct(&mut self, angle: AngularDistribution) {
        self.inner.angle = angle;
    }

    #[setter(angle)]
    pub fn set_angle(&mut self, angle: &PyAny) -> PyResult<()> {
        let typed_angle = angle.extract::<AngularDistribution>()?;
        self.set_angle_direct(typed_angle);
        Ok(())
    }

    pub fn __repr__(&self) -> String {
        format!("IndependentSource(space={:?}, energy={})", 
                self.inner.space, self.inner.energy)
    }
}
