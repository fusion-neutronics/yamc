/// Helper: compute cumulative distribution from PDF (lin-lin, normalized)
fn cumulative_from_pdf(x: &[f64], p: &[f64]) -> Vec<f64> {
    let mut c = Vec::with_capacity(x.len());
    let mut sum = 0.0;
    c.push(0.0);
    for i in 1..x.len() {
        // Trapezoidal rule for lin-lin
        let dx = x[i] - x[i - 1];
        let avg = 0.5 * (p[i] + p[i - 1]);
        sum += avg * dx;
        c.push(sum);
    }
    // Normalize
    if sum > 0.0 {
        for v in &mut c {
            *v /= sum;
        }
    }
    c
}
// Reaction product types and distributions
// Refactored to match OpenMC's file structure:
// - This file contains all type definitions (like reaction_product.h)
// - secondary_*.rs files contain only sampling functions (like secondary_*.cpp)

use serde::{Deserialize, Serialize};
use rand::Rng;

// ============================================================================
// SHARED TYPES (would be in reaction_product.h in OpenMC)
// ============================================================================

/// Tabulated probability distribution with x and p values
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tabulated {
    pub x: Vec<f64>,
    pub p: Vec<f64>,
}

impl Tabulated {
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        if self.x.is_empty() || self.p.is_empty() {
            return 0.0;
        }
        
        let xi = rng.gen::<f64>();
        let mut i = 0;
        for (idx, &prob) in self.p.iter().enumerate() {
            if xi <= prob {
                i = idx;
                break;
            }
        }
        
        if i == 0 {
            return self.x[0];
        }
        if i >= self.x.len() {
            return self.x[self.x.len() - 1];
        }
        
        let f = (xi - self.p[i - 1]) / (self.p[i] - self.p[i - 1]);
        self.x[i - 1] + f * (self.x[i] - self.x[i - 1])
    }
    
    pub fn to_cdf(&self) -> Self {
        if self.p.is_empty() || self.x.is_empty() {
            return self.clone();
        }

        // Use proper trapezoidal integration to convert PDF to CDF
        // This accounts for non-uniform spacing in x values
        let cdf_p = cumulative_from_pdf(&self.x, &self.p);

        Tabulated {
            x: self.x.clone(),
            p: cdf_p,
        }
    }
}

/// 1D tabulated function
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
}

impl Tabulated1D {
    pub fn evaluate(&self, energy: f64) -> f64 {
        match self {
            Tabulated1D::Tabulated1D { x, y, .. } => {
                if x.is_empty() || y.is_empty() {
                    return 0.0;
                }
                
                if energy <= x[0] {
                    return y[0];
                }
                if energy >= x[x.len() - 1] {
                    return y[y.len() - 1];
                }
                
                let i = x.binary_search_by(|&val| val.partial_cmp(&energy).unwrap())
                    .unwrap_or_else(|i| i.saturating_sub(1));
                
                if i >= x.len() - 1 {
                    return y[y.len() - 1];
                }
                
                let f = (energy - x[i]) / (x[i + 1] - x[i]);
                y[i] + f * (y[i + 1] - y[i])
            }
        }
    }
}

/// Tabulated probability distribution for outgoing energy
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum TabulatedProbability {
    #[serde(rename = "Tabulated")]
    Tabulated { x: Vec<f64>, p: Vec<f64> },
}

impl TabulatedProbability {
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        match self {
            TabulatedProbability::Tabulated { x, p } => {
                let tabulated = Tabulated { x: x.clone(), p: p.clone() };
                let cdf = tabulated.to_cdf();
                cdf.sample(rng)
            }
        }
    }
    
    pub fn get_x_values(&self) -> &[f64] {
        match self {
            TabulatedProbability::Tabulated { x, .. } => x,
        }
    }
}

/// Angular distribution
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AngleDistribution {
    pub energy: Vec<f64>,
    pub mu: Vec<Tabulated>,
}

impl AngleDistribution {
    /// Sample scattering cosine from angular distribution
    /// Follows OpenMC's AngleDistribution::sample() approach with stochastic interpolation
    pub fn sample<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> f64 {
        if self.energy.is_empty() || self.mu.is_empty() {
            // No distribution - return isotropic
            return 2.0 * rng.gen::<f64>() - 1.0;
        }

        let n = self.energy.len();

        // Find energy bin and calculate interpolation factor
        // OpenMC: distribution_angle.cpp lines 74-83
        let (i, r) = if incoming_energy < self.energy[0] {
            (0, 0.0)
        } else if incoming_energy > self.energy[n - 1] {
            (n - 2, 1.0)
        } else {
            let idx = self.find_energy_index(incoming_energy);
            let interp = if idx + 1 < n && self.energy[idx + 1] > self.energy[idx] {
                (incoming_energy - self.energy[idx]) / (self.energy[idx + 1] - self.energy[idx])
            } else {
                0.0
            };
            (idx, interp)
        };

        // Stochastic interpolation: sample between the ith and (i+1)th bin
        // OpenMC: distribution_angle.cpp lines 85-87
        let bin = if r > rng.gen::<f64>() { i + 1 } else { i };

        // Ensure bin is within bounds
        let bin = bin.min(self.mu.len() - 1);

        // Sample from the selected distribution
        let mu_sample = self.mu[bin].to_cdf().sample(rng);

        // Make sure mu is in range [-1, 1]
        mu_sample.max(-1.0).min(1.0)
    }

    fn find_energy_index(&self, energy: f64) -> usize {
        // Find lower bound index for interpolation
        for i in 0..self.energy.len() - 1 {
            if energy < self.energy[i + 1] {
                return i;
            }
        }
        self.energy.len() - 2
    }
}

/// Energy distribution types
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum EnergyDistribution {
    LevelInelastic {
        /// Threshold energy for this level
        threshold: f64,
        /// Mass ratio = (A/(A+1))^2 for two-body kinematics
        mass_ratio: f64,
    },
    Tabulated {
        energy: Vec<f64>,
        energy_out: Vec<Vec<f64>>,
    },
    #[serde(rename = "ContinuousTabular")]
    ContinuousTabular {
        energy: Vec<f64>,
        energy_out: Vec<TabulatedProbability>,
    },
}

impl EnergyDistribution {
    pub fn sample<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> f64 {
        match self {
            EnergyDistribution::LevelInelastic { threshold, mass_ratio } => {
                // OpenMC LevelInelastic: E_out = mass_ratio * (E_in - threshold)
                // The threshold is stored such that E_out is valid when E_in > threshold
                // mass_ratio = (A/(A+1))^2 accounts for two-body kinematics
                if incoming_energy > *threshold {
                    mass_ratio * (incoming_energy - threshold)
                } else {
                    // Below threshold - this shouldn't happen in valid physics
                    0.0
                }
            }
            EnergyDistribution::Tabulated { energy, energy_out } => {
                if energy.is_empty() || energy_out.is_empty() {
                    return incoming_energy;
                }
                let i = Self::find_energy_index(incoming_energy, energy);
                if i >= energy_out.len() {
                    return incoming_energy;
                }
                let xs = &energy_out[i];
                if xs.is_empty() {
                    return incoming_energy;
                }
                let idx = rng.gen_range(0..xs.len());
                xs[idx]
            },
            EnergyDistribution::ContinuousTabular { energy, energy_out } => {
                if energy.is_empty() || energy_out.is_empty() {
                    return incoming_energy;
                }

                let n_energy = energy.len();

                // Find energy bin and interpolation factor (OpenMC style)
                let (i, r) = if incoming_energy <= energy[0] {
                    (0, 0.0)
                } else if incoming_energy >= energy[n_energy - 1] {
                    (n_energy.saturating_sub(2), 1.0)
                } else {
                    let idx = Self::find_energy_index(incoming_energy, energy);
                    let interp = if idx + 1 < n_energy && energy[idx + 1] > energy[idx] {
                        (incoming_energy - energy[idx]) / (energy[idx + 1] - energy[idx])
                    } else {
                        0.0
                    };
                    (idx, interp)
                };

                // Check bounds
                if i >= energy_out.len() || i + 1 >= energy_out.len() {
                    if i < energy_out.len() {
                        return energy_out[i].sample(rng);
                    }
                    return incoming_energy;
                }

                // Stochastically choose bin l (either i or i+1)
                let l = if r > rng.gen::<f64>() { i + 1 } else { i };

                // Get E_out bounds for both bins
                let x_i = energy_out[i].get_x_values();
                let x_i1 = energy_out[i + 1].get_x_values();

                if x_i.is_empty() || x_i1.is_empty() {
                    return energy_out[l].sample(rng);
                }

                let e_i_1 = x_i[0];  // min E_out for bin i
                let e_i_k = *x_i.last().unwrap_or(&e_i_1);  // max E_out for bin i
                let e_i1_1 = x_i1[0];  // min E_out for bin i+1
                let e_i1_k = *x_i1.last().unwrap_or(&e_i1_1);  // max E_out for bin i+1

                // Interpolated E_out bounds
                let e_1 = e_i_1 + r * (e_i1_1 - e_i_1);
                let e_k = e_i_k + r * (e_i1_k - e_i_k);

                // Sample from chosen distribution
                let e_out_raw = energy_out[l].sample(rng);

                // Interpolate outgoing energy between incident energy bins (OpenMC style)
                let (e_l_1, e_l_k) = if l == i {
                    (e_i_1, e_i_k)
                } else {
                    (e_i1_1, e_i1_k)
                };

                if e_l_k > e_l_1 && e_k > e_1 {
                    e_1 + (e_out_raw - e_l_1) * (e_k - e_1) / (e_l_k - e_l_1)
                } else {
                    e_out_raw
                }
            }
        }
    }
    
    fn find_energy_index(target_energy: f64, energy_grid: &[f64]) -> usize {
        if target_energy <= energy_grid[0] {
            return 0;
        }
        if target_energy >= energy_grid[energy_grid.len() - 1] {
            return energy_grid.len() - 1;
        }
        energy_grid.binary_search_by(|&val| val.partial_cmp(&target_energy).unwrap())
            .unwrap_or_else(|i| i.saturating_sub(1))
    }
}

// ============================================================================
// MAIN TYPES
// ============================================================================

/// Particle types for reaction products
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum ParticleType {
    #[serde(rename = "neutron")]
    Neutron,
    #[serde(rename = "photon")]
    Photon,
}

/// Angle-energy distribution types
/// Corresponds to OpenMC's AngleEnergy base class with derived types:
/// - UncorrelatedAngleEnergy (secondary_uncorrelated.cpp)
/// - KalbachMann (secondary_kalbach.cpp)
/// - CorrelatedAngleEnergy (secondary_correlated.cpp)
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum AngleEnergyDistribution {
    UncorrelatedAngleEnergy {
        angle: AngleDistribution,
        energy: Option<EnergyDistribution>,
    },
    KalbachMann {
        #[serde(flatten)]
        kalbach: crate::secondary_kalbach::KalbachMann,
    },
    CorrelatedAngleEnergy {
        #[serde(flatten)]
        correlated: crate::secondary_correlated::CorrelatedAngleEnergy,
    },
    /// Evaporation energy distribution (OpenMC-compatible)
    #[serde(rename = "Evaporation")]
    Evaporation {
        theta: Option<Tabulated1D>, // nuclear temperature parameter as function of E
        u: f64, // restriction energy (threshold)
    },
}

impl AngleEnergyDistribution {
    /// Sample outgoing energy and scattering cosine
    /// Delegates to the appropriate secondary distribution sampler
    pub fn sample<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> (f64, f64) {
        match self {
            AngleEnergyDistribution::UncorrelatedAngleEnergy { angle, energy } => {
                crate::secondary_uncorrelated::sample_uncorrelated(incoming_energy, angle, energy, rng)
            },
            AngleEnergyDistribution::KalbachMann { kalbach } => {
                kalbach.sample(incoming_energy, rng)
            },
            AngleEnergyDistribution::CorrelatedAngleEnergy { correlated } => {
                correlated.sample(incoming_energy, rng)
            }
            AngleEnergyDistribution::Evaporation { theta, u } => {
                // OpenMC: sample outgoing energy from evaporation spectrum
                // p(E) ~ exp(-(E-u)/theta(E)), 0 < E_out < E-u
                // Sample: rejection method with two random numbers
                let e_in = incoming_energy;
                let theta_val = if let Some(tab) = theta {
                    tab.evaluate(e_in)
                } else {
                    1.0 // fallback, should not happen
                };
                let y = (e_in - *u) / theta_val;
                let v = 1.0 - (-y).exp();
                let mut e_out = 0.0;
                loop {
                    let xi1 = rng.gen::<f64>();
                    let xi2 = rng.gen::<f64>();
                    let x = -((1.0 - v * xi1) * (1.0 - v * xi2)).ln();
                    if x <= y {
                        e_out = x * theta_val;
                        break;
                    }
                }
                // Isotropic emission
                let mu = 2.0 * rng.gen::<f64>() - 1.0;
                (e_out, mu)
            }
        }
    }
}

/// Yield (multiplicity) of reaction product as function of energy
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum Yield {
    #[serde(rename = "Polynomial")]
    Polynomial { coefficients: Vec<f64> },
    #[serde(rename = "Tabulated")]
    Tabulated { x: Vec<f64>, y: Vec<f64> },
    #[serde(rename = "Tabulated1D")]
    Tabulated1D { 
        x: Vec<f64>, 
        y: Vec<f64>,
        breakpoints: Option<Vec<usize>>,
        interpolation: Option<Vec<u32>>,
    },
}

impl Yield {
    /// Evaluate yield at given energy
    pub fn evaluate(&self, energy: f64) -> f64 {
        match self {
            Yield::Polynomial { coefficients } => {
                if coefficients.is_empty() {
                    return 1.0;
                }
                coefficients.iter().enumerate()
                    .fold(0.0, |acc, (i, &coeff)| acc + coeff * energy.powi(i as i32))
            },
            Yield::Tabulated { x, y } => {
                if x.is_empty() || y.is_empty() {
                    return 1.0;
                }
                if energy <= x[0] {
                    return y[0];
                }
                if energy >= x[x.len() - 1] {
                    return y[y.len() - 1];
                }
                for i in 0..x.len() - 1 {
                    if energy >= x[i] && energy <= x[i + 1] {
                        let f = (energy - x[i]) / (x[i + 1] - x[i]);
                        return y[i] + f * (y[i + 1] - y[i]);
                    }
                }
                1.0
            },
            Yield::Tabulated1D { x, y, .. } => {
                if x.is_empty() || y.is_empty() {
                    return 1.0;
                }
                if energy <= x[0] {
                    return y[0];
                }
                if energy >= x[x.len() - 1] {
                    return y[y.len() - 1];
                }
                for i in 0..x.len() - 1 {
                    if energy >= x[i] && energy <= x[i + 1] {
                        let f = (energy - x[i]) / (x[i + 1] - x[i]);
                        return y[i] + f * (y[i + 1] - y[i]);
                    }
                }
                1.0
            }
        }
    }
}

/// Reaction product with distributions
/// Corresponds to OpenMC's ReactionProduct class
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ReactionProduct {
    pub particle: ParticleType,
    pub emission_mode: String,
    pub decay_rate: f64,
    pub applicability: Vec<Tabulated1D>,
    pub distribution: Vec<AngleEnergyDistribution>,
    #[serde(rename = "yield", default)]
    pub product_yield: Option<Yield>,
}

impl ReactionProduct {
    /// Sample an outgoing particle from this product
    pub fn sample<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> (f64, f64) {
        if self.distribution.is_empty() {
            return (incoming_energy, 2.0 * rng.gen::<f64>() - 1.0);
        }
        
        let distribution_index = if self.distribution.len() == 1 {
            0
        } else {
            self.sample_distribution_index(incoming_energy, rng)
        };
        
        self.distribution[distribution_index].sample(incoming_energy, rng)
    }
    
    fn sample_distribution_index<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> usize {
        if self.applicability.is_empty() || self.applicability.len() != self.distribution.len() {
            return 0;
        }
        
        let mut probs: Vec<f64> = self.applicability.iter()
            .map(|app| app.evaluate(incoming_energy))
            .collect();
        
        let total: f64 = probs.iter().sum();
        if total > 0.0 {
            for p in &mut probs {
                *p /= total;
            }
        } else {
            let n = probs.len() as f64;
            for p in &mut probs {
                *p = 1.0 / n;
            }
        }
        
        let xi = rng.gen::<f64>();
        let mut cumsum = 0.0;
        for (i, &p) in probs.iter().enumerate() {
            cumsum += p;
            if xi <= cumsum {
                return i;
            }
        }
        
        self.distribution.len() - 1
    }
    
    pub fn sample_multiple<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> Vec<(f64, f64)> {
        self.distribution.iter()
            .map(|dist| dist.sample(incoming_energy, rng))
            .collect()
    }
    
    pub fn is_particle_type(&self, particle_type: &ParticleType) -> bool {
        std::mem::discriminant(&self.particle) == std::mem::discriminant(particle_type)
    }
    
    pub fn is_prompt(&self) -> bool {
        self.emission_mode == "prompt"
    }
    
    pub fn is_delayed(&self) -> bool {
        self.emission_mode == "delayed"
    }
    
    pub fn get_decay_rate(&self) -> f64 {
        self.decay_rate
    }
}
