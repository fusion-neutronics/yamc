use crate::reaction_product::{
    ReactionProduct, AngleEnergyDistribution, AngleDistribution, 
    EnergyDistribution, ParticleType, Tabulated
};
use pyo3::prelude::*;
use rand::{thread_rng, Rng};

#[cfg(feature = "pyo3")]
#[pyclass(name = "ReactionProduct")]
#[derive(Clone, Debug)]
pub struct PyReactionProduct {
    pub inner: ReactionProduct,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyReactionProduct {
    #[new]
    fn new(
        particle: String,
        emission_mode: String,
        decay_rate: f64,
    ) -> Self {
        let particle_type = match particle.as_str() {
            "neutron" => ParticleType::Neutron,
            "photon" => ParticleType::Photon,
            _ => ParticleType::Neutron, // Default fallback
        };
        
        Self {
            inner: ReactionProduct {
                particle: particle_type,
                emission_mode,
                decay_rate,
                applicability: vec![],
                distribution: vec![],
                product_yield: None,
            }
        }
    }
    
    /// Sample an outgoing particle from this product
    /// Returns tuple of (outgoing_energy, mu_cosine)
    fn sample(&self, incoming_energy: f64) -> (f64, f64) {
        let mut rng = thread_rng();
        self.inner.sample(incoming_energy, &mut rng)
    }
    
    /// Sample multiple outgoing particles
    /// Returns list of (outgoing_energy, mu_cosine) tuples
    fn sample_multiple(&self, incoming_energy: f64) -> Vec<(f64, f64)> {
        let mut rng = thread_rng();
        self.inner.sample_multiple(incoming_energy, &mut rng)
    }
    
    /// Check if this product represents a specific particle type
    fn is_particle_type(&self, particle_name: &str) -> bool {
        let particle_type = match particle_name {
            "neutron" => ParticleType::Neutron,
            "photon" => ParticleType::Photon,
            _ => return false,
        };
        self.inner.is_particle_type(&particle_type)
    }
    
    /// Check if this is a prompt emission
    fn is_prompt(&self) -> bool {
        self.inner.is_prompt()
    }
    
    /// Check if this is a delayed emission
    fn is_delayed(&self) -> bool {
        self.inner.is_delayed()
    }
    
    /// Get the decay rate for delayed neutron precursors
    fn get_decay_rate(&self) -> f64 {
        self.inner.get_decay_rate()
    }
    
    /// Get particle type as string
    #[getter]
    fn particle(&self) -> String {
        match self.inner.particle {
            ParticleType::Neutron => "neutron".to_string(),
            ParticleType::Photon => "photon".to_string(),
        }
    }
    
    /// Get emission mode
    #[getter]
    fn emission_mode(&self) -> String {
        self.inner.emission_mode.clone()
    }
    
    /// Get number of distributions
    #[getter]
    fn num_distributions(&self) -> usize {
        self.inner.distribution.len()
    }
    
    /// Get distribution types as strings
    #[getter]
    fn distribution_types(&self) -> Vec<String> {
        self.inner.distribution.iter().map(|dist| {
            match dist {
                crate::reaction_product::AngleEnergyDistribution::UncorrelatedAngleEnergy { .. } => "UncorrelatedAngleEnergy".to_string(),
                crate::reaction_product::AngleEnergyDistribution::KalbachMann { .. } => "KalbachMann".to_string(),
                crate::reaction_product::AngleEnergyDistribution::CorrelatedAngleEnergy { .. } => "CorrelatedAngleEnergy".to_string(),
            }
        }).collect()
    }
}

#[cfg(feature = "pyo3")]
impl PyReactionProduct {
    /// Create PyReactionProduct from a Rust ReactionProduct
    pub fn from_reaction_product(product: ReactionProduct, _py: pyo3::Python) -> pyo3::PyResult<Self> {
        Ok(PyReactionProduct { inner: product })
    }
}

#[cfg(feature = "pyo3")]
#[pyclass]
#[derive(Clone)]
pub struct PyAngleDistribution {
    pub inner: AngleDistribution,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyAngleDistribution {
    /// Sample scattering cosine (mu) for a given incoming energy
    fn sample(&self, incoming_energy: f64) -> f64 {
        let mut rng = thread_rng();
        self.inner.sample(incoming_energy, &mut rng)
    }
    
    /// Get energy grid
    #[getter]
    fn energy(&self) -> Vec<f64> {
        self.inner.energy.clone()
    }
    
    /// Get number of energy points
    #[getter]
    fn num_energy_points(&self) -> usize {
        self.inner.energy.len()
    }
}

#[cfg(feature = "pyo3")]
#[pyclass]
#[derive(Clone)]
pub struct PyTabulated {
    pub inner: Tabulated,
}

#[cfg(feature = "pyo3")]
#[pymethods]
impl PyTabulated {
    #[new]
    fn new(x: Vec<f64>, p: Vec<f64>) -> Self {
        Self {
            inner: Tabulated { x, p }
        }
    }
    
    /// Sample from the tabulated distribution
    fn sample(&self) -> f64 {
        let mut rng = thread_rng();
        self.inner.sample(&mut rng)
    }
    
    /// Convert to CDF and return as new PyTabulated
    fn to_cdf(&self) -> PyTabulated {
        PyTabulated {
            inner: self.inner.to_cdf()
        }
    }
    
    /// Get x values
    #[getter]
    fn x(&self) -> Vec<f64> {
        self.inner.x.clone()
    }
    
    /// Get p values
    #[getter]
    fn p(&self) -> Vec<f64> {
        self.inner.p.clone()
    }
}

/// Convenience function to sample isotropic scattering
#[cfg(feature = "pyo3")]
#[pyfunction]
pub fn sample_isotropic() -> f64 {
    let mut rng = thread_rng();
    2.0 * rng.gen::<f64>() - 1.0
}

/// Convenience function to create tabulated distribution and sample from it
#[cfg(feature = "pyo3")]
#[pyfunction]
pub fn sample_tabulated(x: Vec<f64>, p: Vec<f64>) -> f64 {
    let tabulated = Tabulated { x, p };
    let cdf = tabulated.to_cdf();
    let mut rng = thread_rng();
    cdf.sample(&mut rng)
}

/// Test function to create a simple reaction product for testing
#[cfg(feature = "pyo3")]
#[pyfunction]
pub fn create_test_reaction_product() -> PyReactionProduct {
    // Create a simple test product with elastic scattering (energy unchanged)
    
    // Create simple angular distribution (isotropic-ish)
    let mu_dist = Tabulated {
        x: vec![-1.0, 0.0, 1.0],
        p: vec![0.33, 0.67, 1.0], // CDF values
    };
    
    let angle_dist = AngleDistribution {
        energy: vec![1e5, 1e6, 1e7], // 0.1, 1, 10 MeV
        mu: vec![mu_dist.clone(), mu_dist.clone(), mu_dist.clone()],
    };
    
    // Create uncorrelated angle-energy distribution
    let angle_energy_dist = AngleEnergyDistribution::UncorrelatedAngleEnergy {
        angle: angle_dist,
        energy: None, // No energy change (elastic)
    };
    
    let product = ReactionProduct {
        particle: ParticleType::Neutron,
        emission_mode: "prompt".to_string(),
        decay_rate: 0.0,
        applicability: vec![],
        distribution: vec![angle_energy_dist],
        product_yield: None,
    };
    
    PyReactionProduct { inner: product }
}