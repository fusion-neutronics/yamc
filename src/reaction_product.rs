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
pub struct AngleDistribution {
    pub energy: Vec<f64>,
    pub mu: Vec<Tabulated>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnergyDistribution {
    pub energy: Vec<f64>,
    pub energy_out: Vec<Vec<f64>>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum AngleEnergyDistribution {
    UncorrelatedAngleEnergy {
        angle: AngleDistribution,
        energy: Option<EnergyDistribution>,
    },
    // Add other distribution types as needed (CorrelatedAngleEnergy, KalbachMann, etc.)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionProduct {
    /// Type of particle (e.g., neutron, photon)
    pub particle: ParticleType,
    /// Emission mode (prompt, delayed, total)
    pub emission_mode: String,
    /// Decay rate (for delayed neutron precursors) in [1/s]
    pub decay_rate: f64,
    /// Applicability of distribution (empty for now)
    pub applicability: Vec<serde_json::Value>, // Placeholder for applicability data
    /// Distributions of energy and angle of product
    pub distribution: Vec<AngleEnergyDistribution>,
}
