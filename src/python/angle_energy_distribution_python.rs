use pyo3::prelude::*;
use crate::reaction_product::{AngleDistribution, EnergyDistribution, Tabulated1D, TabulatedProbability};

#[pyclass(name = "UncorrelatedAngleEnergy")]
pub struct PyUncorrelatedAngleEnergy {
    #[pyo3(get)]
    pub angle: PyAngleDistribution,
    #[pyo3(get)]
    pub energy: Option<PyEnergyDistribution>,
}

#[pymethods]
impl PyUncorrelatedAngleEnergy {
    #[new]
    fn new(angle: PyAngleDistribution, energy: Option<PyEnergyDistribution>) -> Self {
        PyUncorrelatedAngleEnergy { angle, energy }
    }

    fn __repr__(&self) -> String {
        format!("UncorrelatedAngleEnergy(angle={:?}, energy={:?})", 
                self.angle.energy.len(), self.energy.is_some())
    }
}

#[pyclass(name = "KalbachMann")]
pub struct PyKalbachMann {
    #[pyo3(get)]
    pub energy: Vec<f64>,
    #[pyo3(get)]
    pub energy_out: Vec<pyo3::PyObject>,
    #[pyo3(get)]
    pub slope: Vec<pyo3::PyObject>,
}

#[pymethods]
impl PyKalbachMann {
    #[new]
    fn new(energy: Vec<f64>, energy_out: Vec<pyo3::PyObject>, slope: Vec<pyo3::PyObject>) -> Self {
        PyKalbachMann { energy, energy_out, slope }
    }

    fn __repr__(&self) -> String {
        format!("KalbachMann(energy=[{} points], energy_out=[{} points], slope=[{} points])", 
                self.energy.len(), self.energy_out.len(), self.slope.len())
    }
}

// Existing wrappers for AngleDistribution and EnergyDistribution
#[pyclass(name = "AngleDistribution")]
#[derive(Clone)]
pub struct PyAngleDistribution {
    #[pyo3(get)]
    pub energy: Vec<f64>,
    #[pyo3(get)]
    pub mu: Vec<PyTabulated>,
}

#[pyclass(name = "EnergyDistribution")]
#[derive(Clone)]
pub struct PyEnergyDistribution {
    pub distribution_type: String,
    pub data: pyo3::PyObject,
}

#[pyclass(name = "Tabulated")]
#[derive(Clone)]
pub struct PyTabulated {
    #[pyo3(get)]
    pub x: Vec<f64>,
    #[pyo3(get)]
    pub p: Vec<f64>,
}

#[pymethods]
impl PyTabulated {
    #[new]
    fn new(x: Vec<f64>, p: Vec<f64>) -> Self {
        PyTabulated { x, p }
    }

    fn __repr__(&self) -> String {
        format!("Tabulated(x=[{} points], p=[{} points])", self.x.len(), self.p.len())
    }
}

impl From<crate::reaction_product::Tabulated> for PyTabulated {
    fn from(tab: crate::reaction_product::Tabulated) -> Self {
        PyTabulated { x: tab.x, p: tab.p }
    }
}

#[pymethods]
impl PyAngleDistribution {
    #[new]
    fn new(energy: Vec<f64>, mu: Vec<PyTabulated>) -> Self {
        PyAngleDistribution { energy, mu }
    }

    fn __repr__(&self) -> String {
        format!(
            "AngleDistribution(energy=[{} points], mu=[{} distributions])",
            self.energy.len(),
            self.mu.len()
        )
    }
}

impl From<AngleDistribution> for PyAngleDistribution {
    fn from(angle_dist: AngleDistribution) -> Self {
        PyAngleDistribution {
            energy: angle_dist.energy,
            mu: angle_dist.mu.into_iter().map(PyTabulated::from).collect(),
        }
    }
}

#[pymethods]
impl PyEnergyDistribution {
    #[new]
    fn new(distribution_type: String, data: PyObject) -> Self {
        PyEnergyDistribution { distribution_type, data }
    }

    #[getter]
    fn distribution_type(&self) -> &str {
        &self.distribution_type
    }

    #[getter]
    fn data(&self) -> PyObject {
        self.data.clone()
    }

    fn __repr__(&self) -> String {
        format!("EnergyDistribution(type='{}')", self.distribution_type)
    }
}

impl PyEnergyDistribution {
    pub fn from_energy_distribution(energy_dist: EnergyDistribution, py: Python) -> PyResult<Self> {
        use pyo3::types::PyDict;

        let (dist_type, data) = match energy_dist {
            EnergyDistribution::LevelInelastic {} => {
                let dict = PyDict::new(py);
                ("LevelInelastic".to_string(), dict.into())
            }
            EnergyDistribution::Tabulated { energy, energy_out } => {
                let dict = PyDict::new(py);
                dict.set_item("energy", energy)?;
                dict.set_item("energy_out", energy_out)?;
                ("Tabulated".to_string(), dict.into())
            }
        };

        Ok(PyEnergyDistribution {
            distribution_type: dist_type,
            data,
        })
    }
}

// Python wrappers for Tabulated1D
#[pyclass(name = "Tabulated1D")]
#[derive(Clone)]
pub struct PyTabulated1D {
    #[pyo3(get)]
    pub x: Vec<f64>,
    #[pyo3(get)]
    pub y: Vec<f64>,
    #[pyo3(get)]
    pub breakpoints: Vec<i32>,
    #[pyo3(get)]
    pub interpolation: Vec<i32>,
}

#[pymethods]
impl PyTabulated1D {
    #[new]
    fn new(x: Vec<f64>, y: Vec<f64>, breakpoints: Vec<i32>, interpolation: Vec<i32>) -> Self {
        PyTabulated1D { x, y, breakpoints, interpolation }
    }

    fn __repr__(&self) -> String {
        format!(
            "Tabulated1D(x=[{} points], y=[{} points], breakpoints=[{}], interpolation=[{}])",
            self.x.len(),
            self.y.len(),
            self.breakpoints.len(),
            self.interpolation.len()
        )
    }
}

impl From<Tabulated1D> for PyTabulated1D {
    fn from(tab: Tabulated1D) -> Self {
        match tab {
            Tabulated1D::Tabulated1D { x, y, breakpoints, interpolation } => {
                PyTabulated1D { x, y, breakpoints, interpolation }
            }
        }
    }
}

// Add From<TabulatedProbability> for PyTabulated to handle energy_out conversions
impl From<TabulatedProbability> for PyTabulated {
    fn from(tab: TabulatedProbability) -> Self {
        match tab {
            TabulatedProbability::Tabulated { x, p } => {
                PyTabulated { x, p }
            }
        }
    }
}
