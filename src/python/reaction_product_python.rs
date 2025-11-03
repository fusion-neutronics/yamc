#[cfg(feature = "pyo3")]
use crate::python::angle_energy_distribution_python::{PyAngleDistribution, PyEnergyDistribution, PyTabulated};
use serde_json::Value;
#[cfg(feature = "pyo3")]
use crate::python::serde_json_utils_python::*;
use crate::reaction_product::{
    AngleDistribution, AngleEnergyDistribution, EnergyDistribution, ParticleType, ReactionProduct,
    Tabulated,
};
#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use pyo3::types::PyDict;

#[cfg(feature = "pyo3")]
#[pyclass(name = "ReactionProduct")]
#[derive(Clone, Debug)]
pub struct PyReactionProduct {
    #[pyo3(get)]
    pub particle: String,
    #[pyo3(get)]
    pub emission_mode: Option<String>,
    #[pyo3(get)]
    pub decay_rate: Option<f64>,
    #[pyo3(get)]
    pub applicability: Vec<PyObject>,
    #[pyo3(get)]
    pub distribution: Vec<PyObject>,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyReactionProduct {
    #[new]
    fn new(
        particle: String,
        emission_mode: Option<String>,
        decay_rate: Option<f64>,
        applicability: Vec<PyObject>,
        distribution: Vec<PyObject>,
    ) -> Self {
        PyReactionProduct {
            particle,
            emission_mode,
            decay_rate,
            applicability,
            distribution,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "ReactionProduct(particle='{}', distribution=[{} items])",
            self.particle,
            self.distribution.len()
        )
    }
}

#[cfg(feature = "pyo3")]



#[cfg(feature = "pyo3")]
#[pyclass(name = "AngleEnergyDistribution")]
#[derive(Clone, Debug)]
pub struct PyAngleEnergyDistribution {
    pub distribution_type: String,
    pub data: PyObject,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyAngleEnergyDistribution {
    #[getter]
    fn distribution_type(&self) -> &str {
        &self.distribution_type
    }

    #[getter]
    fn data(&self) -> PyObject {
        self.data.clone()
    }

    fn __repr__(&self) -> String {
        format!("AngleEnergyDistribution(type='{}')", self.distribution_type)
    }
}

#[cfg(feature = "pyo3")]
impl PyAngleEnergyDistribution {
    fn from_angle_energy_distribution(
        angle_energy_dist: AngleEnergyDistribution,
        py: Python,
    ) -> PyResult<Self> {
        use pyo3::types::PyDict;

        let (dist_type, data) = match angle_energy_dist {
            AngleEnergyDistribution::UncorrelatedAngleEnergy { angle, energy } => {
                let dict = PyDict::new(py);
                dict.set_item("angle", Py::new(py, PyAngleDistribution::from(angle))?)?;
                if let Some(energy_dist) = energy {
                    dict.set_item("energy", Py::new(py, PyEnergyDistribution::from_energy_distribution(energy_dist, py)?)?)?;
                }
                ("UncorrelatedAngleEnergy".to_string(), dict.into())
            }
            AngleEnergyDistribution::KalbachMann {
                energy,
                energy_out,
                slope,
            } => {
                let dict = PyDict::new(py);
                dict.set_item("energy", energy)?;
                dict.set_item("energy_out", energy_out.into_iter().map(|v| value_to_pyobject(&v, py)).collect::<Vec<_>>())?;
                dict.set_item("slope", slope.into_iter().map(|v| value_to_pyobject(&v, py)).collect::<Vec<_>>())?;
                ("KalbachMann".to_string(), dict.into())
            }
        };
        Ok(PyAngleEnergyDistribution {
            distribution_type: dist_type,
            data,
        })
    }
}

#[cfg(feature = "pyo3")]
impl PyReactionProduct {
    pub fn from_reaction_product(product: ReactionProduct, py: Python) -> PyResult<Self> {
        let particle = match product.particle {
            ParticleType::Neutron => "neutron".to_string(),
            ParticleType::Photon => "photon".to_string(),
        };

        let mut distribution = Vec::new();
        for dist in product.distribution {
            match dist {
                AngleEnergyDistribution::UncorrelatedAngleEnergy { angle, energy } => {
                    let py_angle = PyAngleDistribution {
                        energy: angle.energy.clone(),
                        mu: angle.mu.into_iter().map(|tab| PyTabulated { x: tab.x, p: tab.p }).collect(),
                    };
                    let py_energy = match energy {
                        Some(e) => Some(PyEnergyDistribution::from_energy_distribution(e, py)?),
                        None => None,
                    };
                    let py_dist = pyo3::Py::new(py, crate::python::angle_energy_distribution_python::PyUncorrelatedAngleEnergy {
                        angle: py_angle,
                        energy: py_energy,
                    })?;
                    distribution.push(py_dist.into_py(py));
                }
                AngleEnergyDistribution::KalbachMann { energy, energy_out, slope } => {
                    let py_dist = pyo3::Py::new(py, crate::python::angle_energy_distribution_python::PyKalbachMann {
                        energy,
                        energy_out: energy_out.into_iter().map(|v| value_to_pyobject(&v, py)).collect(),
                        slope: slope.into_iter().map(|v| value_to_pyobject(&v, py)).collect(),
                    })?;
                    distribution.push(py_dist.into_py(py));
                }
            }
        }

        // Convert applicability (currently just raw JSON values)
        let applicability = product
            .applicability
            .into_iter()
                .map(|val| value_to_pyobject(&val, py))
            .collect();

        Ok(PyReactionProduct {
            particle,
            emission_mode: Some(product.emission_mode),
            decay_rate: Some(product.decay_rate),
            applicability,
            distribution,
        })
    }
}