// Uncorrelated angle-energy distributions
// Corresponds to OpenMC's secondary_uncorrelated.cpp

// use serde::{Deserialize, Serialize};
use rand::Rng;
use crate::reaction_product::{EnergyDistribution, AngleDistribution};

/// Sample from uncorrelated angle-energy distribution
/// This function corresponds to OpenMC's UncorrelatedAngleEnergy::sample()
///
/// The angle and energy are sampled independently:
/// - If angle distribution exists, sample mu from it
/// - Otherwise use isotropic scattering
/// - If energy distribution exists, sample E_out from it
/// - Otherwise return None for energy (caller must handle kinematics)
pub fn sample_uncorrelated<R: Rng>(
    incoming_energy: f64,
    angle: &AngleDistribution,
    energy: &Option<EnergyDistribution>,
    rng: &mut R,
) -> (f64, f64) {
    // Sample cosine of scattering angle
    let mu = angle.sample(incoming_energy, rng);

    // Sample outgoing energy if distribution exists
    // If no energy distribution, return incoming energy as placeholder
    // (caller should apply two-body kinematics for elastic scattering)
    let e_out = energy.as_ref()
        .map(|e| e.sample(incoming_energy, rng))
        .unwrap_or(incoming_energy);

    (e_out, mu)
}

/// Sample from uncorrelated angle-energy for elastic scattering
/// Applies two-body kinematics to compute outgoing energy from angle
pub fn sample_elastic_uncorrelated<R: Rng>(
    incoming_energy: f64,
    angle: &AngleDistribution,
    awr: f64,  // atomic weight ratio (target mass / neutron mass)
    rng: &mut R,
) -> (f64, f64) {
    // Sample cosine of scattering angle in center-of-mass frame
    let mu_cm = angle.sample(incoming_energy, rng);

    // Two-body elastic kinematics (target at rest approximation)
    // In CM frame: E_out/E_in = (1 + alpha + (1 - alpha)*mu_cm) / 2
    // where alpha = ((A-1)/(A+1))^2
    let alpha = ((awr - 1.0) / (awr + 1.0)).powi(2);
    let e_out = incoming_energy * (1.0 + alpha + (1.0 - alpha) * mu_cm) / 2.0;

    // Convert mu from CM to lab frame
    // mu_lab = (1 + A*mu_cm) / sqrt(1 + A^2 + 2*A*mu_cm)
    let denom = (1.0 + awr * awr + 2.0 * awr * mu_cm).sqrt();
    let mu_lab = (1.0 + awr * mu_cm) / denom;

    (e_out, mu_lab)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    use crate::reaction_product::Tabulated;
    
    #[test]
    fn test_isotropic_angle_fallback() {
        let angle = AngleDistribution {
            energy: vec![],
            mu: vec![],
        };
        
        let mut rng = StdRng::seed_from_u64(42);
        let mu = angle.sample(1e6, &mut rng);
        
        // Should be isotropic: -1 <= mu <= 1
        assert!(mu >= -1.0 && mu <= 1.0);
    }
    
    #[test]
    fn test_uncorrelated_elastic() {
        let angle = AngleDistribution {
            energy: vec![1e5, 1e6],
            mu: vec![
                Tabulated { x: vec![-1.0, 0.0, 1.0], p: vec![0.0, 0.5, 1.0] },
                Tabulated { x: vec![-1.0, 0.0, 1.0], p: vec![0.0, 0.5, 1.0] },
            ],
        };
        
        let mut rng = StdRng::seed_from_u64(42);
        let (e_out, mu) = sample_uncorrelated(5e5, &angle, &None, &mut rng);
        
        // No energy distribution, so elastic
        assert_eq!(e_out, 5e5);
        // Mu should be sampled from distribution
        assert!(mu >= -1.0 && mu <= 1.0);
    }
}
