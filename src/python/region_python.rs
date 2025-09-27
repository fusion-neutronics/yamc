impl PyRegionExpr {
    pub fn to_region_expr(&self) -> crate::region::RegionExpr {
        match self {
            PyRegionExpr::Halfspace(hs) => {
                // Convert PyHalfspace to HalfspaceType
                Python::with_gil(|py| {
                    let surface = hs.surface.as_ref(py);
                    let rust_surface = surface.borrow().inner.clone();
                    let arc_surface = std::sync::Arc::new(rust_surface);
                    if hs.is_above {
                        crate::region::RegionExpr::Halfspace(crate::region::HalfspaceType::Above(
                            arc_surface,
                        ))
                    } else {
                        crate::region::RegionExpr::Halfspace(crate::region::HalfspaceType::Below(
                            arc_surface,
                        ))
                    }
                })
            }
            PyRegionExpr::Union(a, b) => crate::region::RegionExpr::Union(
                Box::new(a.to_region_expr()),
                Box::new(b.to_region_expr()),
            ),
            PyRegionExpr::Intersection(a, b) => crate::region::RegionExpr::Intersection(
                Box::new(a.to_region_expr()),
                Box::new(b.to_region_expr()),
            ),
            PyRegionExpr::Complement(inner) => {
                crate::region::RegionExpr::Complement(Box::new(inner.to_region_expr()))
            }
        }
    }
}
use pyo3::prelude::*;
// ...existing code...

// ...existing code...
use crate::surface_python::PySurface;

#[pyclass(name = "Region")]
#[derive(Clone)]
pub struct PyRegion {
    pub expr: PyRegionExpr,
    pub region: crate::region::Region,
}

#[derive(Clone)]
pub enum PyRegionExpr {
    Halfspace(PyHalfspace),
    Union(Box<PyRegionExpr>, Box<PyRegionExpr>),
    Intersection(Box<PyRegionExpr>, Box<PyRegionExpr>),
    Complement(Box<PyRegionExpr>),
}

#[pymethods]
impl PyRegion {
    fn __invert__(self_: PyRef<'_, Self>) -> PyResult<Self> {
        let expr = PyRegionExpr::Complement(Box::new(self_.expr.clone()));
        let region = crate::region::Region {
            expr: expr.to_region_expr(),
        };
        Ok(PyRegion { expr, region })
    }

    pub fn contains(&self, point: (f64, f64, f64)) -> bool {
        self.region.contains(point)
    }

    pub fn bounding_box(&self) -> PyBoundingBox {
        use crate::bounding_box::BoundingBox;
        let bbox: BoundingBox = self.region.bounding_box();
        PyBoundingBox {
            lower_left: bbox.lower_left,
            upper_right: bbox.upper_right,
            center: [
                (bbox.lower_left[0] + bbox.upper_right[0]) / 2.0,
                (bbox.lower_left[1] + bbox.upper_right[1]) / 2.0,
                (bbox.lower_left[2] + bbox.upper_right[2]) / 2.0,
            ],
            width: [
                bbox.upper_right[0] - bbox.lower_left[0],
                bbox.upper_right[1] - bbox.lower_left[1],
                bbox.upper_right[2] - bbox.lower_left[2],
            ],
        }
    }

    fn __and__(&self, other: &PyAny) -> PyResult<PyRegion> {
        let expr = if let Ok(other_region) = other.extract::<PyRef<PyRegion>>() {
            PyRegionExpr::Intersection(
                Box::new(self.expr.clone()),
                Box::new(other_region.expr.clone()),
            )
        } else if let Ok(other_halfspace) = other.extract::<PyRef<PyHalfspace>>() {
            PyRegionExpr::Intersection(
                Box::new(self.expr.clone()),
                Box::new(PyRegionExpr::Halfspace(other_halfspace.clone())),
            )
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Operand must be PyRegion or PyHalfspace",
            ));
        };
        let region = crate::region::Region {
            expr: expr.to_region_expr(),
        };
        Ok(PyRegion { expr, region })
    }

    fn __or__(&self, other: &PyAny) -> PyResult<PyRegion> {
        let expr = if let Ok(other_region) = other.extract::<PyRef<PyRegion>>() {
            PyRegionExpr::Union(
                Box::new(self.expr.clone()),
                Box::new(other_region.expr.clone()),
            )
        } else if let Ok(other_halfspace) = other.extract::<PyRef<PyHalfspace>>() {
            PyRegionExpr::Union(
                Box::new(self.expr.clone()),
                Box::new(PyRegionExpr::Halfspace(other_halfspace.clone())),
            )
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Operand must be PyRegion or PyHalfspace",
            ));
        };
        let region = crate::region::Region {
            expr: expr.to_region_expr(),
        };
        Ok(PyRegion { expr, region })
    }
}

#[pyclass]
pub struct PyBoundingBox {
    #[pyo3(get)]
    pub lower_left: [f64; 3],
    #[pyo3(get)]
    pub upper_right: [f64; 3],
    #[pyo3(get)]
    pub center: [f64; 3],
    #[pyo3(get)]
    pub width: [f64; 3],
}

#[pymethods]
impl PyBoundingBox {
    // Optionally, add __repr__ for pretty printing
    fn __repr__(&self) -> String {
        format!(
            "BoundingBox(lower_left={:?}, upper_right={:?}, center={:?}, width={:?})",
            self.lower_left, self.upper_right, self.center, self.width
        )
    }
}

#[pyclass]
#[derive(Clone)]
pub struct PyHalfspace {
    pub surface: Py<PySurface>,
    pub is_above: bool,
}

#[pymethods]
impl PyHalfspace {
    #[staticmethod]
    pub fn new_above(surface: Py<PySurface>) -> Self {
        PyHalfspace {
            surface,
            is_above: true,
        }
    }
    #[staticmethod]
    pub fn new_below(surface: Py<PySurface>) -> Self {
        PyHalfspace {
            surface,
            is_above: false,
        }
    }
    fn __neg__(slf: PyRef<'_, Self>) -> PyResult<Self> {
        Ok(PyHalfspace {
            surface: slf.surface.clone(),
            is_above: false,
        })
    }
    fn __pos__(slf: PyRef<'_, Self>) -> PyResult<Self> {
        Ok(PyHalfspace {
            surface: slf.surface.clone(),
            is_above: true,
        })
    }
    fn __invert__(slf: PyRef<'_, Self>) -> PyResult<PyRegion> {
        let expr = PyRegionExpr::Complement(Box::new(PyRegionExpr::Halfspace(slf.clone())));
        let region = crate::region::Region {
            expr: expr.to_region_expr(),
        };
        Ok(PyRegion { expr, region })
    }
    pub fn contains(&self, point: (f64, f64, f64)) -> bool {
        Python::with_gil(|py| {
            let surface = self.surface.as_ref(py);
            if self.is_above {
                surface.borrow().evaluate(point) > 0.0
            } else {
                surface.borrow().evaluate(point) < 0.0
            }
        })
    }
    pub fn bounding_box(&self) -> PyBoundingBox {
        Python::with_gil(|py| {
            let surface = self.surface.as_ref(py);
            match &surface.borrow().inner.kind {
                crate::surface::SurfaceKind::Plane { a, b, c, d } => {
                    let mut lower = [f64::NEG_INFINITY; 3];
                    let mut upper = [f64::INFINITY; 3];
                    if *a == 1.0 && *b == 0.0 && *c == 0.0 {
                        if self.is_above {
                            lower[0] = *d;
                        } else {
                            upper[0] = *d;
                        }
                    } else if *a == 0.0 && *b == 1.0 && *c == 0.0 {
                        if self.is_above {
                            lower[1] = *d;
                        } else {
                            upper[1] = *d;
                        }
                    } else if *a == 0.0 && *b == 0.0 && *c == 1.0 {
                        if self.is_above {
                            lower[2] = *d;
                        } else {
                            upper[2] = *d;
                        }
                    }
                    PyBoundingBox {
                        lower_left: lower,
                        upper_right: upper,
                        center: [0.0, 0.0, 0.0],
                        width: [0.0, 0.0, 0.0],
                    }
                }
                crate::surface::SurfaceKind::Sphere { x0, y0, z0, radius } => PyBoundingBox {
                    lower_left: [*x0 - *radius, *y0 - *radius, *z0 - *radius],
                    upper_right: [*x0 + *radius, *y0 + *radius, *z0 + *radius],
                    center: [*x0, *y0, *z0],
                    width: [2.0 * *radius, 2.0 * *radius, 2.0 * *radius],
                },
                _ => PyBoundingBox {
                    lower_left: [f64::NEG_INFINITY; 3],
                    upper_right: [f64::INFINITY; 3],
                    center: [0.0, 0.0, 0.0],
                    width: [0.0, 0.0, 0.0],
                },
            }
        })
    }
    fn __and__(&self, other: &PyAny) -> PyResult<PyRegion> {
        let expr = if let Ok(other_halfspace) = other.extract::<PyRef<PyHalfspace>>() {
            PyRegionExpr::Intersection(
                Box::new(PyRegionExpr::Halfspace(self.clone())),
                Box::new(PyRegionExpr::Halfspace(other_halfspace.clone())),
            )
        } else if let Ok(other_region) = other.extract::<PyRef<PyRegion>>() {
            PyRegionExpr::Intersection(
                Box::new(PyRegionExpr::Halfspace(self.clone())),
                Box::new(other_region.expr.clone()),
            )
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Operand must be PyRegion or PyHalfspace",
            ));
        };
        let region = crate::region::Region {
            expr: expr.to_region_expr(),
        };
        Ok(PyRegion { expr, region })
    }
    fn __or__(&self, other: &PyAny) -> PyResult<PyRegion> {
        let expr = if let Ok(other_halfspace) = other.extract::<PyRef<PyHalfspace>>() {
            PyRegionExpr::Union(
                Box::new(PyRegionExpr::Halfspace(self.clone())),
                Box::new(PyRegionExpr::Halfspace(other_halfspace.clone())),
            )
        } else if let Ok(other_region) = other.extract::<PyRef<PyRegion>>() {
            PyRegionExpr::Union(
                Box::new(PyRegionExpr::Halfspace(self.clone())),
                Box::new(other_region.expr.clone()),
            )
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Operand must be PyRegion or PyHalfspace",
            ));
        };
        let region = crate::region::Region {
            expr: expr.to_region_expr(),
        };
        Ok(PyRegion { expr, region })
    }
}

impl PyRegionExpr {
    pub fn evaluate_contains(&self, point: (f64, f64, f64)) -> bool {
        match self {
            PyRegionExpr::Halfspace(hs) => Python::with_gil(|py| {
                let surface = hs.surface.as_ref(py);
                if hs.is_above {
                    surface.borrow().evaluate(point) > 0.0
                } else {
                    surface.borrow().evaluate(point) < 0.0
                }
            }),
            PyRegionExpr::Union(a, b) => a.evaluate_contains(point) || b.evaluate_contains(point),
            PyRegionExpr::Intersection(a, b) => {
                a.evaluate_contains(point) && b.evaluate_contains(point)
            }
            PyRegionExpr::Complement(inner) => !inner.evaluate_contains(point),
        }
    }
    // ...existing code...
}
