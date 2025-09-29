
use crate::cell::Cell;
use crate::python::region_python::PyRegion;
use crate::material::Material;
use pyo3::prelude::*;
use pyo3::types::PyAny;

// Use the PyMaterial definition from material_python.rs

#[pyclass(name = "Cell")]
#[derive(Clone)]
pub struct PyCell {
    pub inner: Cell,
}

#[pymethods]
impl PyCell {
    /// Compute the distance to the closest surface from a point in a direction (if any)
    pub fn distance_to_surface(&self, point: (f64, f64, f64), direction: (f64, f64, f64)) -> Option<f64> {
        let point_arr = [point.0, point.1, point.2];
        let dir_arr = [direction.0, direction.1, direction.2];
        self.inner.distance_to_surface(point_arr, dir_arr)
    }

    /// Return the closest surface (as PySurface) and distance, or None if no intersection
    pub fn closest_surface(&self, point: (f64, f64, f64), direction: (f64, f64, f64)) -> Option<(crate::python::surface_python::PySurface, f64)> {
        let point_arr = [point.0, point.1, point.2];
        let dir_arr = [direction.0, direction.1, direction.2];
        if let Some(surf_arc) = self.inner.closest_surface(point_arr, dir_arr) {
            if let Some(dist) = surf_arc.distance_to_surface(point_arr, dir_arr) {
                return Some((crate::python::surface_python::PySurface { inner: (*surf_arc).clone() }, dist));
            }
        }
        None
    }
    #[new]
    #[pyo3(signature = (region, cell_id=0, name=None, fill=None))]
    pub fn new(region: PyRegion, cell_id: u32, name: Option<String>, fill: Option<crate::python::material_python::PyMaterial>) -> Self {
        let material = fill.map(|mat| std::sync::Arc::new(std::sync::Mutex::new(mat.internal.clone())));
        PyCell {
            inner: Cell::new(cell_id, region.region, name, material),
        }
    }

    #[getter]
    pub fn cell_id(&self) -> u32 {
        self.inner.cell_id
    }

    #[getter]
    pub fn name(&self) -> Option<String> {
        self.inner.name.clone()
    }
    #[getter]
    pub fn fill(&self) -> Option<crate::python::material_python::PyMaterial> {
        self.inner.material.as_ref().map(|arc| {
            let mat = arc.lock().unwrap().clone();
            crate::python::material_python::PyMaterial { internal: mat }
        })
    }


    pub fn contains(&self, x: f64, y: f64, z: f64) -> bool {
        self.inner.contains((x, y, z))
    }
}
