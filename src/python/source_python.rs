use crate::source::IndependentSource;
use pyo3::prelude::*;
use pyo3::types::PyAny;

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
        angle: Option<&PyAny>
    ) -> PyResult<Self> {
        let mut instance = Self { inner: IndependentSource::new() };
        
        if let Some(val) = space { instance.set_space(val); }
        if let Some(val) = energy { instance.set_energy(val); }
        if let Some(val) = angle { instance.set_angle(val)?; }
        
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
    pub fn angle(&self) -> PyObject {
        Python::with_gil(|py| self.convert_angle_to_python(py))
    }
    
    fn convert_angle_to_python(&self, py: Python) -> PyObject {
        use crate::python::stats_python::{PyIsotropic, PyMonodirectional};
        use crate::stats::{Isotropic, Monodirectional};
        use std::any::Any;
        
        let angle_ref = self.inner.angle.as_ref();
        
        match (
            angle_ref.as_any().downcast_ref::<Isotropic>(),
            angle_ref.as_any().downcast_ref::<Monodirectional>()
        ) {
            (Some(iso), _) => PyIsotropic { inner: iso.clone() }.into_py(py),
            (_, Some(mono)) => PyMonodirectional { inner: mono.clone() }.into_py(py),
            _ => format!("AngularDistribution(type={})", angle_ref.type_name()).to_object(py)
        }
    }

    #[setter(angle)]
    pub fn set_angle(&mut self, angle: &PyAny) -> PyResult<()> {
        use crate::python::stats_python::{PyIsotropic, PyMonodirectional};
        
        // Try to extract as Monodirectional first
        if let Ok(mono) = angle.extract::<PyMonodirectional>() {
            self.inner.angle = Box::new(mono.inner.clone());
        } 
        // Try to extract as Isotropic
        else if let Ok(iso) = angle.extract::<PyIsotropic>() {
            self.inner.angle = Box::new(iso.inner.clone());
        }
        // If neither, return error
        else {
            return Err(PyErr::new::<pyo3::exceptions::PyTypeError, _>(
                "angle must be an Isotropic or Monodirectional distribution"
            ));
        }
        Ok(())
    }

    pub fn __repr__(&self) -> String {
        format!("IndependentSource(space={:?}, energy={})", 
                self.inner.space, self.inner.energy)
    }
}
