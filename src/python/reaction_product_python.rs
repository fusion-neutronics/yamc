use crate::reaction_product::{
    AngleDistribution, AngleEnergyDistribution, EnergyDistribution, ParticleType, ReactionProduct,
    Tabulated,
};
#[cfg(feature = "pyo3")]
use pyo3::prelude::*;
#[cfg(feature = "pyo3")]
use pyo3::types::PyDict;

#[cfg(feature = "pyo3")]
#[pyclass(name = "Tabulated")]
#[derive(Clone, Debug)]
pub struct PyTabulated {
    #[pyo3(get)]
    pub x: Vec<f64>,
    #[pyo3(get)]
    pub p: Vec<f64>,
}

#[cfg(feature = "pyo3")]
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

#[cfg(feature = "pyo3")]
impl From<Tabulated> for PyTabulated {
    fn from(tab: Tabulated) -> Self {
        PyTabulated { x: tab.x, p: tab.p }
    }
}

#[cfg(feature = "pyo3")]
#[pyclass(name = "AngleDistribution")]
#[derive(Clone, Debug)]
pub struct PyAngleDistribution {
    #[pyo3(get)]
    pub energy: Vec<f64>,
    #[pyo3(get)]
    pub mu: Vec<PyTabulated>,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyAngleDistribution {
    fn __repr__(&self) -> String {
        format!(
            "AngleDistribution(energy=[{} points], mu=[{} distributions])",
            self.energy.len(),
            self.mu.len()
        )
    }
}

#[cfg(feature = "pyo3")]
impl From<AngleDistribution> for PyAngleDistribution {
    fn from(angle_dist: AngleDistribution) -> Self {
        PyAngleDistribution {
            energy: angle_dist.energy,
            mu: angle_dist.mu.into_iter().map(PyTabulated::from).collect(),
        }
    }
}

#[cfg(feature = "pyo3")]
#[pyclass(name = "EnergyDistribution")]
#[derive(Clone, Debug)]
pub struct PyEnergyDistribution {
    pub distribution_type: String,
    pub data: PyObject,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyEnergyDistribution {
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

#[cfg(feature = "pyo3")]
impl PyEnergyDistribution {
    fn from_energy_distribution(energy_dist: EnergyDistribution, py: Python) -> PyResult<Self> {
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
                dict.set_item("angle", PyAngleDistribution::from(angle))?;
                if let Some(energy_dist) = energy {
                    dict.set_item("energy", PyEnergyDistribution::from_energy_distribution(energy_dist, py)?)?;
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
                dict.set_item("energy_out", energy_out)?;
                dict.set_item("slope", slope)?;
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
#[pyclass(name = "ReactionProduct")]
#[derive(Clone, Debug)]
pub struct PyReactionProduct {
    #[pyo3(get)]
    pub particle: String,
    #[pyo3(get)]
    pub emission_mode: String,
    #[pyo3(get)]
    pub decay_rate: f64,
    #[pyo3(get)]
    pub applicability: Vec<PyObject>,
    #[pyo3(get)]
    pub distribution: Vec<PyAngleEnergyDistribution>,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyReactionProduct {
    fn __repr__(&self) -> String {
        format!(
            "ReactionProduct(particle='{}', emission_mode='{}', distributions={})",
            self.particle,
            self.emission_mode,
            self.distribution.len()
        )
    }
}

#[cfg(feature = "pyo3")]
impl PyReactionProduct {
    pub fn from_reaction_product(product: ReactionProduct, py: Python) -> PyResult<Self> {
        let particle = match product.particle {
            ParticleType::Neutron => "neutron".to_string(),
            ParticleType::Photon => "photon".to_string(),
        };

        let distribution = product
            .distribution
            .into_iter()
            .map(|dist| PyAngleEnergyDistribution::from_angle_energy_distribution(dist, py))
            .collect::<PyResult<Vec<_>>>()?;

        // Convert applicability (currently just raw JSON values)
        let applicability = product
            .applicability
            .into_iter()
            .map(|val| val.into_py(py))
            .collect();

        Ok(PyReactionProduct {
            particle,
            emission_mode: product.emission_mode,
            decay_rate: product.decay_rate,
            applicability,
            distribution,
        })
    }
}