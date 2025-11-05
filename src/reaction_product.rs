use serde::{Deserialize, Serialize};
use rand::Rng;

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

impl Tabulated {
    /// Sample from a tabulated probability distribution using inverse CDF method
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        if self.x.is_empty() || self.p.is_empty() {
            return 0.0;
        }
        
        let xi = rng.gen::<f64>();
        
        // Find the interval in the CDF
        let mut i = 0;
        for (idx, &prob) in self.p.iter().enumerate() {
            if xi <= prob {
                i = idx;
                break;
            }
        }
        
        // Handle edge cases
        if i == 0 {
            return self.x[0];
        }
        if i >= self.x.len() {
            return self.x[self.x.len() - 1];
        }
        
        // Linear interpolation between points
        let f = (xi - self.p[i - 1]) / (self.p[i] - self.p[i - 1]);
        self.x[i - 1] + f * (self.x[i] - self.x[i - 1])
    }
    
    /// Convert probability density to cumulative distribution
    pub fn to_cdf(&self) -> Self {
        if self.p.is_empty() {
            return self.clone();
        }
        
        let mut cdf_p = Vec::with_capacity(self.p.len());
        let mut cumsum = 0.0;
        
        for &prob in &self.p {
            cumsum += prob;
            cdf_p.push(cumsum);
        }
        
        // Normalize to ensure CDF goes to 1.0
        let total = cdf_p[cdf_p.len() - 1];
        if total > 0.0 {
            for p in &mut cdf_p {
                *p /= total;
            }
        }
        
        Tabulated {
            x: self.x.clone(),
            p: cdf_p,
        }
    }
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

impl Tabulated1D {
    /// Evaluate the tabulated function at a given energy
    pub fn evaluate(&self, energy: f64) -> f64 {
        match self {
            Tabulated1D::Tabulated1D { x, y, .. } => {
                if x.is_empty() || y.is_empty() {
                    return 0.0;
                }
                
                // Find the energy bracket
                if energy <= x[0] {
                    return y[0];
                }
                if energy >= x[x.len() - 1] {
                    return y[y.len() - 1];
                }
                
                // Binary search for the interval
                let i = x.binary_search_by(|&val| val.partial_cmp(&energy).unwrap())
                    .unwrap_or_else(|i| i.saturating_sub(1));
                
                if i >= x.len() - 1 {
                    return y[y.len() - 1];
                }
                
                // Linear interpolation
                let f = (energy - x[i]) / (x[i + 1] - x[i]);
                y[i] + f * (y[i + 1] - y[i])
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum TabulatedProbability {
    #[serde(rename = "Tabulated")]
    Tabulated { x: Vec<f64>, p: Vec<f64> },
    // Add other probability distribution types as needed
}

impl TabulatedProbability {
    /// Sample from the tabulated probability distribution
    pub fn sample<R: Rng>(&self, rng: &mut R) -> f64 {
        match self {
            TabulatedProbability::Tabulated { x, p } => {
                let tabulated = Tabulated { x: x.clone(), p: p.clone() };
                let cdf = tabulated.to_cdf();
                cdf.sample(rng)
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AngleDistribution {
    pub energy: Vec<f64>,
    pub mu: Vec<Tabulated>,
}

impl AngleDistribution {
    /// Sample scattering cosine (mu) for a given incoming energy
    pub fn sample<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> f64 {
        if self.energy.is_empty() || self.mu.is_empty() {
            // Isotropic scattering as fallback
            return 2.0 * rng.gen::<f64>() - 1.0;
        }
        
        // Find the energy bracket
        let i = self.find_energy_index(incoming_energy);
        
        if i >= self.mu.len() {
            let mu_sample = self.mu[self.mu.len() - 1].to_cdf().sample(rng);
            return mu_sample.max(-1.0).min(1.0);
        }
        
        // If we're at a tabulated energy point, sample directly
        if i < self.energy.len() && (incoming_energy - self.energy[i]).abs() < 1e-10 {
            let mu_sample = self.mu[i].to_cdf().sample(rng);
            return mu_sample.max(-1.0).min(1.0);
        }
        
        // Interpolation between energy points
        if i > 0 && i < self.mu.len() {
            // For simplicity, just use the closest distribution
            // A more sophisticated approach would interpolate the CDFs
            let mu_sample = if i < self.energy.len() {
                let f = (incoming_energy - self.energy[i - 1]) / (self.energy[i] - self.energy[i - 1]);
                if f < 0.5 {
                    self.mu[i - 1].to_cdf().sample(rng)
                } else {
                    self.mu[i].to_cdf().sample(rng)
                }
            } else {
                self.mu[i - 1].to_cdf().sample(rng)
            };
            
            // Ensure mu is within valid bounds
            mu_sample.max(-1.0).min(1.0)
        } else {
            // Use the closest available distribution and clamp result
            let mu_sample = self.mu[i.min(self.mu.len() - 1)].to_cdf().sample(rng);
            mu_sample.max(-1.0).min(1.0)
        }
    }
    
    fn find_energy_index(&self, energy: f64) -> usize {
        if energy <= self.energy[0] {
            return 0;
        }
        if energy >= self.energy[self.energy.len() - 1] {
            return self.energy.len() - 1;
        }
        
        // Binary search for energy bracket
        self.energy.binary_search_by(|&val| val.partial_cmp(&energy).unwrap())
            .unwrap_or_else(|i| i)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "type")]
pub enum EnergyDistribution {
    LevelInelastic {
        // Level inelastic scattering - outgoing energy is deterministic
    },
    Tabulated {
        energy: Vec<f64>,
        energy_out: Vec<Vec<f64>>,
    },
    // Add other energy distribution types as needed
}

impl EnergyDistribution {
    /// Sample outgoing energy for a given incoming energy
    pub fn sample<R: Rng>(&self, incoming_energy: f64, _rng: &mut R) -> f64 {
        match self {
            EnergyDistribution::LevelInelastic { .. } => {
                // For level inelastic, outgoing energy equals incoming energy
                // minus the Q-value (would need Q-value data for proper implementation)
                incoming_energy
            },
            EnergyDistribution::Tabulated { energy, energy_out } => {
                if energy.is_empty() || energy_out.is_empty() {
                    return incoming_energy;
                }
                
                // Find energy bracket
                let i = self.find_energy_index(incoming_energy, energy);
                
                if i >= energy_out.len() {
                    return incoming_energy;
                }
                
                // For now, return the first energy in the distribution
                // A complete implementation would sample from the energy_out distribution
                if !energy_out[i].is_empty() {
                    energy_out[i][0]
                } else {
                    incoming_energy
                }
            }
        }
    }
    
    fn find_energy_index(&self, target_energy: f64, energy_grid: &[f64]) -> usize {
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
    CorrelatedAngleEnergy {
        energy: Vec<f64>,
        energy_out: Vec<Vec<f64>>,
        mu: Vec<Vec<TabulatedProbability>>,
        // Note: Following OpenMC's CorrelatedAngleEnergy format where mu contains
        // the angular distributions as Tabulated objects, not separate distribution field
    },
    // Add other distribution types as needed
}

impl AngleEnergyDistribution {
    /// Sample outgoing energy and scattering cosine
    pub fn sample<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> (f64, f64) {
        match self {
            AngleEnergyDistribution::UncorrelatedAngleEnergy { angle, energy } => {
                // Sample angle and energy independently
                let mu = angle.sample(incoming_energy, rng);
                let e_out = energy.as_ref()
                    .map(|e| e.sample(incoming_energy, rng))
                    .unwrap_or(incoming_energy);
                (e_out, mu)
            },
            AngleEnergyDistribution::KalbachMann { energy, energy_out, slope } => {
                // Kalbach-Mann correlated angle-energy sampling
                self.sample_kalbach_mann(incoming_energy, energy, energy_out, slope, rng)
            },
            AngleEnergyDistribution::CorrelatedAngleEnergy { energy, energy_out, mu } => {
                // Correlated angle-energy sampling following OpenMC's implementation
                self.sample_correlated_angle_energy_openmc(incoming_energy, energy, energy_out, mu, rng)
            }
        }
    }
    
    fn sample_kalbach_mann<R: Rng>(
        &self,
        incoming_energy: f64,
        energy_grid: &[f64],
        energy_out: &[TabulatedProbability],
        slope: &[Tabulated1D],
        rng: &mut R,
    ) -> (f64, f64) {
        // Find energy bracket
        let i = self.find_energy_index(incoming_energy, energy_grid);
        
        if i >= energy_out.len() || i >= slope.len() {
            // Fallback to isotropic
            return (incoming_energy, 2.0 * rng.gen::<f64>() - 1.0);
        }
        
        // Sample outgoing energy from distribution
        let e_out = energy_out[i].sample(rng);
        
        // Sample angle using Kalbach-Mann systematics
        let a = slope[i].evaluate(e_out);
        
        // Kalbach-Mann angular distribution: f(mu) ~ exp(a*mu)
        // Use rejection sampling for simplicity
        let mu = if a.abs() < 1e-6 {
            // Isotropic case
            2.0 * rng.gen::<f64>() - 1.0
        } else {
            // Sample from exp(a*mu) on [-1, 1]
            let r = rng.gen::<f64>();
            if a > 0.0 {
                ((1.0 - r) * (-a).exp() + r * a.exp()).ln() / a
            } else {
                -((1.0 - r) * a.exp() + r * (-a).exp()).ln() / a
            }
        };
        
        // Ensure mu is in [-1, 1]
        let mu = mu.max(-1.0).min(1.0);
        
        (e_out, mu)
    }
    
    fn sample_correlated_angle_energy_openmc<R: Rng>(
        &self,
        incoming_energy: f64,
        energy_grid: &[f64],
        energy_out_grid: &[Vec<f64>], 
        mu_distributions: &[Vec<TabulatedProbability>],
        rng: &mut R,
    ) -> (f64, f64) {
        // Following OpenMC's CorrelatedAngleEnergy implementation
        // Find incoming energy bracket
        let energy_index = self.find_energy_index(incoming_energy, energy_grid);
        
        if energy_index >= energy_out_grid.len() || energy_index >= mu_distributions.len() {
            // Outside tabulated range - fallback to isotropic elastic
            return (incoming_energy, 2.0 * rng.gen::<f64>() - 1.0);
        }
        
        let e_out_grid = &energy_out_grid[energy_index];
        let angular_distributions = &mu_distributions[energy_index];
        
        if e_out_grid.is_empty() || angular_distributions.is_empty() {
            // No data available - fallback to isotropic elastic
            return (incoming_energy, 2.0 * rng.gen::<f64>() - 1.0);
        }
        
        // Sample outgoing energy index (uniform for simplicity - should be from energy distribution)
        let e_out_index = if e_out_grid.len() == 1 {
            0
        } else {
            rng.gen_range(0..e_out_grid.len().min(angular_distributions.len()))
        };
        
        let e_out = e_out_grid[e_out_index];
        
        // Sample scattering cosine from the corresponding angular distribution
        let mu = if e_out_index < angular_distributions.len() {
            angular_distributions[e_out_index].sample(rng)
        } else {
            // Fallback to isotropic
            2.0 * rng.gen::<f64>() - 1.0
        };
        
        // Ensure mu is in valid range [-1, 1]
        let mu = mu.max(-1.0).min(1.0);
        
        (e_out, mu)
    }
    
    fn find_energy_index(&self, target_energy: f64, energy_grid: &[f64]) -> usize {
        if energy_grid.is_empty() {
            return 0;
        }
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

impl ReactionProduct {
    /// Sample an outgoing particle from this product
    /// Returns (outgoing_energy, mu_cosine)
    pub fn sample<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> (f64, f64) {
        if self.distribution.is_empty() {
            // No distribution data, assume elastic scattering
            return (incoming_energy, 2.0 * rng.gen::<f64>() - 1.0);
        }
        
        // If multiple distributions, sample which one to use based on applicability
        let distribution_index = if self.distribution.len() == 1 {
            0
        } else {
            self.sample_distribution_index(incoming_energy, rng)
        };
        
        let distribution = &self.distribution[distribution_index];
        distribution.sample(incoming_energy, rng)
    }
    
    /// Sample which distribution to use based on applicability probabilities
    fn sample_distribution_index<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> usize {
        if self.applicability.is_empty() || self.applicability.len() != self.distribution.len() {
            // No applicability data or mismatch, choose first distribution
            return 0;
        }
        
        // Evaluate applicability probabilities at incoming energy
        let mut probs: Vec<f64> = self.applicability.iter()
            .map(|app| app.evaluate(incoming_energy))
            .collect();
        
        // Normalize probabilities
        let total: f64 = probs.iter().sum();
        if total > 0.0 {
            for p in &mut probs {
                *p /= total;
            }
        } else {
            // Equal probability fallback
            let n = probs.len() as f64;
            for p in &mut probs {
                *p = 1.0 / n;
            }
        }
        
        // Sample from discrete distribution
        let xi = rng.gen::<f64>();
        let mut cumsum = 0.0;
        for (i, &p) in probs.iter().enumerate() {
            cumsum += p;
            if xi <= cumsum {
                return i;
            }
        }
        
        // Fallback to last distribution
        self.distribution.len() - 1
    }
    
    /// Sample multiple outgoing particles if this product represents multiple emissions
    pub fn sample_multiple<R: Rng>(&self, incoming_energy: f64, rng: &mut R) -> Vec<(f64, f64)> {
        let mut results = Vec::new();
        
        for distribution in &self.distribution {
            let (e_out, mu) = distribution.sample(incoming_energy, rng);
            results.push((e_out, mu));
        }
        
        results
    }
    
    /// Check if this product represents a specific particle type
    pub fn is_particle_type(&self, particle_type: &ParticleType) -> bool {
        std::mem::discriminant(&self.particle) == std::mem::discriminant(particle_type)
    }
    
    /// Check if this is a prompt emission (vs delayed)
    pub fn is_prompt(&self) -> bool {
        self.emission_mode == "prompt"
    }
    
    /// Check if this is a delayed emission
    pub fn is_delayed(&self) -> bool {
        self.emission_mode == "delayed"
    }
    
    /// Get the decay rate for delayed neutron precursors
    pub fn get_decay_rate(&self) -> f64 {
        self.decay_rate
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    
    #[test]
    fn test_tabulated_cdf_conversion() {
        // Test PDF to CDF conversion
        let pdf = Tabulated {
            x: vec![-1.0, -0.5, 0.0, 0.5, 1.0],
            p: vec![0.1, 0.2, 0.4, 0.2, 0.1], // PDF values
        };
        
        let cdf = pdf.to_cdf();
        
        // Check that CDF is properly normalized
        assert!((cdf.p.last().unwrap() - 1.0).abs() < 1e-10);
        
        // Check CDF is monotonically increasing
        for i in 1..cdf.p.len() {
            assert!(cdf.p[i] >= cdf.p[i-1]);
        }
    }
    
    #[test]
    fn test_tabulated_sampling() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Simple uniform-ish distribution
        let cdf = Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.33, 0.67, 1.0],
        };
        
        // Sample multiple values
        let samples: Vec<f64> = (0..1000).map(|_| cdf.sample(&mut rng)).collect();
        
        // Check samples are within bounds
        for sample in &samples {
            assert!(*sample >= -1.0 && *sample <= 1.0);
        }
        
        // Check rough distribution (this is probabilistic)
        let mean: f64 = samples.iter().sum::<f64>() / samples.len() as f64;
        assert!(mean.abs() < 0.5); // Should be roughly centered around 0 (relaxed tolerance)
    }
    
    #[test]
    fn test_tabulated1d_evaluation() {
        let tab1d = Tabulated1D::Tabulated1D {
            x: vec![1e5, 1e6, 1e7],
            y: vec![10.0, 5.0, 1.0],
            breakpoints: vec![],
            interpolation: vec![],
        };
        
        // Test exact points
        assert!((tab1d.evaluate(1e5) - 10.0).abs() < 1e-10);
        assert!((tab1d.evaluate(1e6) - 5.0).abs() < 1e-10);
        assert!((tab1d.evaluate(1e7) - 1.0).abs() < 1e-10);
        
        // Test interpolation
        let mid_val = tab1d.evaluate(5.5e5); // Halfway between 1e5 and 1e6
        assert!(mid_val > 5.0 && mid_val < 10.0);
        
        // Test extrapolation (should return boundary values)
        assert!((tab1d.evaluate(1e4) - 10.0).abs() < 1e-10); // Below range
        assert!((tab1d.evaluate(1e8) - 1.0).abs() < 1e-10);  // Above range
    }
    
    #[test]
    fn test_angle_distribution_sampling() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Create simple angular distribution
        let mu_dist = Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.2, 0.6, 1.0], // Slightly forward-peaked
        };
        
        let angle_dist = AngleDistribution {
            energy: vec![1e6, 1e7],
            mu: vec![mu_dist.clone(), mu_dist.clone()],
        };
        
        // Test sampling at different energies
        let mu1 = angle_dist.sample(1e6, &mut rng);
        let mu2 = angle_dist.sample(5e6, &mut rng); // Interpolated
        let mu3 = angle_dist.sample(1e7, &mut rng);
        
        // All should be valid cosines
        assert!(mu1 >= -1.0 && mu1 <= 1.0);
        assert!(mu2 >= -1.0 && mu2 <= 1.0);
        assert!(mu3 >= -1.0 && mu3 <= 1.0);
    }
    
    #[test]
    fn test_energy_distribution_level_inelastic() {
        let mut rng = StdRng::seed_from_u64(42);
        
        let energy_dist = EnergyDistribution::LevelInelastic {};
        
        let incoming_energy = 14e6; // 14 MeV
        let outgoing_energy = energy_dist.sample(incoming_energy, &mut rng);
        
        // For level inelastic without Q-value, energy should be conserved
        assert!((outgoing_energy - incoming_energy).abs() < 1e-10);
    }
    
    #[test]
    fn test_uncorrelated_angle_energy_distribution() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Create angular distribution
        let mu_dist = Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.25, 0.5, 1.0],
        };
        
        let angle_dist = AngleDistribution {
            energy: vec![1e6],
            mu: vec![mu_dist],
        };
        
        // Create angle-energy distribution (elastic - no energy change)
        let dist = AngleEnergyDistribution::UncorrelatedAngleEnergy {
            angle: angle_dist,
            energy: None,
        };
        
        let incoming_energy = 14e6;
        let (e_out, mu) = dist.sample(incoming_energy, &mut rng);
        
        // Energy should be unchanged (elastic)
        assert!((e_out - incoming_energy).abs() < 1e-10);
        
        // Mu should be valid cosine
        assert!(mu >= -1.0 && mu <= 1.0);
    }
    
    #[test]
    fn test_kalbach_mann_distribution() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Create Kalbach-Mann distribution
        let energy_out_dist = TabulatedProbability::Tabulated {
            x: vec![0.1e6, 1e6, 10e6],
            p: vec![0.3, 0.7, 1.0],
        };
        
        let slope = Tabulated1D::Tabulated1D {
            x: vec![0.1e6, 1e6, 10e6],
            y: vec![2.0, 1.0, 0.5], // Kalbach-Mann slope parameter
            breakpoints: vec![],
            interpolation: vec![],
        };
        
        let dist = AngleEnergyDistribution::KalbachMann {
            energy: vec![14e6],
            energy_out: vec![energy_out_dist],
            slope: vec![slope],
        };
        
        let (e_out, mu) = dist.sample(14e6, &mut rng);
        
        // Energy should be sampled from distribution (not necessarily = incoming)
        assert!(e_out > 0.0);
        
        // Mu should be valid cosine
        assert!(mu >= -1.0 && mu <= 1.0);
    }
    
    #[test]
    fn test_reaction_product_sampling() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Create a complete reaction product
        let mu_dist = Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.33, 0.67, 1.0],
        };
        
        let angle_dist = AngleDistribution {
            energy: vec![1e6, 1e7],
            mu: vec![mu_dist.clone(), mu_dist.clone()],
        };
        
        let angle_energy_dist = AngleEnergyDistribution::UncorrelatedAngleEnergy {
            angle: angle_dist,
            energy: Some(EnergyDistribution::LevelInelastic {}),
        };
        
        let product = ReactionProduct {
            particle: ParticleType::Neutron,
            emission_mode: "prompt".to_string(),
            decay_rate: 0.0,
            applicability: vec![],
            distribution: vec![angle_energy_dist],
        };
        
        let incoming_energy = 14e6;
        let (e_out, mu) = product.sample(incoming_energy, &mut rng);
        
        // Check outputs are reasonable
        assert!(e_out > 0.0);
        assert!(mu >= -1.0 && mu <= 1.0);
        
        // Test utility methods
        assert!(product.is_particle_type(&ParticleType::Neutron));
        assert!(!product.is_particle_type(&ParticleType::Photon));
        assert!(product.is_prompt());
        assert!(!product.is_delayed());
        assert!((product.get_decay_rate() - 0.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_reaction_product_multiple_distributions() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Create applicability function (energy-dependent probability)
        let applicability1 = Tabulated1D::Tabulated1D {
            x: vec![1e5, 1e6, 1e7],
            y: vec![1.0, 0.5, 0.1], // Decreases with energy
            breakpoints: vec![],
            interpolation: vec![],
        };
        
        let applicability2 = Tabulated1D::Tabulated1D {
            x: vec![1e5, 1e6, 1e7],
            y: vec![0.1, 0.5, 1.0], // Increases with energy
            breakpoints: vec![],
            interpolation: vec![],
        };
        
        // Create two different distributions
        let mu_dist1 = Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.5, 0.75, 1.0], // Backward-peaked
        };
        
        let mu_dist2 = Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.25, 0.5, 1.0], // Forward-peaked
        };
        
        let angle_dist1 = AngleDistribution {
            energy: vec![1e6],
            mu: vec![mu_dist1],
        };
        
        let angle_dist2 = AngleDistribution {
            energy: vec![1e6],
            mu: vec![mu_dist2],
        };
        
        let dist1 = AngleEnergyDistribution::UncorrelatedAngleEnergy {
            angle: angle_dist1,
            energy: None,
        };
        
        let dist2 = AngleEnergyDistribution::UncorrelatedAngleEnergy {
            angle: angle_dist2,
            energy: None,
        };
        
        let product = ReactionProduct {
            particle: ParticleType::Neutron,
            emission_mode: "prompt".to_string(),
            decay_rate: 0.0,
            applicability: vec![applicability1, applicability2],
            distribution: vec![dist1, dist2],
        };
        
        // Sample at different energies
        let (e_out_low, mu_low) = product.sample(1e5, &mut rng);   // Should favor dist1
        let (e_out_high, mu_high) = product.sample(1e7, &mut rng); // Should favor dist2
        
        // Both should produce valid results
        assert!(e_out_low > 0.0 && mu_low >= -1.0 && mu_low <= 1.0);
        assert!(e_out_high > 0.0 && mu_high >= -1.0 && mu_high <= 1.0);
    }
    
    #[test]
    fn test_correlated_angle_energy_distribution() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Create correlated angle-energy distribution
        let energy_out_dist1 = TabulatedProbability::Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.2, 0.6, 1.0],
        };
        
        let energy_out_dist2 = TabulatedProbability::Tabulated {
            x: vec![-1.0, 0.5, 1.0],
            p: vec![0.1, 0.4, 1.0],
        };
        
        // Create angular distribution as TabulatedProbability
        let mu_dist = TabulatedProbability::Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.3, 0.6, 1.0],  // Must be cumulative
        };
        
        let dist = AngleEnergyDistribution::CorrelatedAngleEnergy {
            energy: vec![1e6],
            energy_out: vec![vec![0.5e6, 1.5e6]],
            mu: vec![vec![mu_dist]],
        };
        
        let (e_out, mu) = dist.sample(1e6, &mut rng);
        
        // Energy should be from the tabulated values
        assert!(e_out > 0.0);
        
        // Mu should be valid cosine
        assert!(mu >= -1.0 && mu <= 1.0);
    }
    
    #[test]
    fn test_isotropic_fallback() {
        let mut rng = StdRng::seed_from_u64(42);
        
        // Empty product should fall back to isotropic scattering
        let product = ReactionProduct {
            particle: ParticleType::Neutron,
            emission_mode: "prompt".to_string(),
            decay_rate: 0.0,
            applicability: vec![],
            distribution: vec![],
        };
        
        let (e_out, mu) = product.sample(14e6, &mut rng);
        
        // Energy should be conserved (elastic fallback)
        assert!((e_out - 14e6).abs() < 1e-10);
        
        // Mu should be isotropic (between -1 and 1)
        assert!(mu >= -1.0 && mu <= 1.0);
    }
}
