use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ParticleType {
    #[serde(rename = "neutron")]
    Neutron,
    #[serde(rename = "photon")]
    Photon,
    // Add more types as needed
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tabulated {
    pub x: Vec<f64>,
    pub p: Vec<f64>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum Tabulated1D {
    #[serde(rename = "Tabulated1D")]
    Tabulated1D {
        x: Vec<f64>,
        y: Vec<f64>,
        breakpoints: Vec<i32>,
        interpolation: Vec<i32>,
    },
    // Add other 1D distribution types as needed
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum TabulatedProbability {
    #[serde(rename = "Tabulated")]
    Tabulated { x: Vec<f64>, p: Vec<f64> },
    // Add other probability distribution types as needed
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AngleDistribution {
    pub energy: Vec<f64>,
    pub mu: Vec<Tabulated>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum EnergyDistribution {
    LevelInelastic {
        // Level inelastic scattering - no additional data needed
    },
    Tabulated {
        energy: Vec<f64>,
        energy_out: Vec<Vec<f64>>,
    },
    // Add other energy distribution types as needed
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum AngleEnergyDistribution {
    UncorrelatedAngleEnergy {
        angle: AngleDistribution,
        energy: Option<EnergyDistribution>,
    },
    KalbachMann {
        energy: Vec<f64>,
        energy_out: Vec<TabulatedProbability>,
        slope: Vec<Tabulated1D>,
    },
    // Add other distribution types as needed (CorrelatedAngleEnergy, etc.)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionProduct {
    /// Type of particle (e.g., neutron, photon)
    pub particle: ParticleType,
    /// Emission mode (prompt, delayed, total)
    pub emission_mode: String,
    /// Decay rate (for delayed neutron precursors) in [1/s]
    pub decay_rate: f64,
    /// Applicability of distribution (energy-dependent probability)
    pub applicability: Vec<Tabulated1D>,
    /// Distributions of energy and angle of product
    pub distribution: Vec<AngleEnergyDistribution>,
}
