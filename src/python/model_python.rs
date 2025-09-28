use crate::model::Model;
use crate::python::geometry_python::PyGeometry;
use crate::python::materials_python::PyMaterials;
use crate::python::settings_python::PySettings;
use crate::python::source_python::PySource;
use pyo3::prelude::*;

#[pyclass(name = "Model")]
#[derive(Clone)]
pub struct PyModel {
    pub inner: Model,
}

#[pymethods]
impl PyModel {
    #[new]
    pub fn new(geometry: PyGeometry, materials: PyMaterials, settings: PySettings) -> Self {
        PyModel {
            inner: Model {
                geometry: geometry.inner.clone(),
                materials: materials.inner.clone(),
                settings: settings.inner.clone(),
            },
        }
    }
    #[getter]
    pub fn geometry(&self) -> PyGeometry {
        PyGeometry {
            inner: self.inner.geometry.clone(),
        }
    }
    #[getter]
    pub fn materials(&self) -> PyMaterials {
        PyMaterials {
            inner: self.inner.materials.clone(),
        }
    }
    #[getter]
    pub fn settings(&self) -> PySettings {
        PySettings {
            inner: self.inner.settings.clone(),
        }
    }

    pub fn run(&self) {
        self.inner.run();
    }
}
