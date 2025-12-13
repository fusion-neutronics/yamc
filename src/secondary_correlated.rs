// Correlated angle-energy distribution sampling
// Corresponds to OpenMC's secondary_correlated.cpp

use serde::{Deserialize, Serialize};
use rand::Rng;
use crate::reaction_product::TabulatedProbability;

/// Sample from correlated angle-energy distribution
/// 
/// This function corresponds to OpenMC's CorrelatedAngleEnergy::sample()
/// 
/// In a correlated distribution, the outgoing angle depends on the outgoing energy.
/// The data structure stores:
/// 1. A distribution of outgoing energies for each incoming energy
/// 2. For each outgoing energy, a distribution of scattering angles (mu)
/// 
/// The sampling procedure:
/// 1. Find the bracket for incoming energy
/// 2. Sample outgoing energy from the distribution
/// 3. Find the bracket for the sampled outgoing energy
/// 4. Sample scattering angle from the corresponding mu distribution
/// 5. Interpolate if needed
/// 
/// # Arguments
/// * `incoming_energy` - Incident particle energy (eV)
/// * `energy_grid` - Tabulated incident energy points
/// * `energy_out_distributions` - Outgoing energy probability distributions at each incident energy
/// * `mu_distributions` - Scattering angle distributions for each (E_in, E_out) combination
/// * `rng` - Random number generator
/// 
/// # Returns
/// Tuple of (outgoing_energy, mu_cosine)
pub fn sample_correlated_angle_energy<R: Rng>(
    incoming_energy: f64,
    energy_grid: &[f64],
    energy_out_distributions: &[TabulatedProbability],
    mu_distributions: &[Vec<TabulatedProbability>],
    rng: &mut R,
) -> (f64, f64) {
    // Find incoming energy bracket
    let energy_index = find_energy_index(incoming_energy, energy_grid);
    
    if energy_index >= energy_out_distributions.len() || energy_index >= mu_distributions.len() {
        // Outside tabulated range - fallback to isotropic elastic
        return (incoming_energy, 2.0 * rng.gen::<f64>() - 1.0);
    }
    
    // Sample outgoing energy from the tabulated distribution at this incoming energy
    let e_out = energy_out_distributions[energy_index].sample(rng);
    
    // Get the angular distributions corresponding to this incoming energy
    let angular_distributions = &mu_distributions[energy_index];
    
    if angular_distributions.is_empty() {
        // No angular data available - fallback to isotropic
        return (e_out, 2.0 * rng.gen::<f64>() - 1.0);
    }
    
    // Get outgoing energy grid from energy_out distribution
    let e_out_grid = energy_out_distributions[energy_index].get_x_values();
    
    // Number of angular distributions available
    let n_ang = angular_distributions.len();
    
    if n_ang == 1 {
        // Only one angular distribution available
        let mu = angular_distributions[0].sample(rng);
        return (e_out, mu.max(-1.0).min(1.0));
    }
    
    // The angular distributions correspond to the outgoing energy bins in e_out_grid
    // Find the bracket for the sampled outgoing energy
    let mut idx = 0;
    if e_out_grid.len() >= 2 && n_ang >= 2 {
        // Find the energy bracket
        for i in 0..e_out_grid.len().saturating_sub(1).min(n_ang - 1) {
            if e_out >= e_out_grid[i] && e_out <= e_out_grid[i + 1] {
                idx = i;
                break;
            }
        }
        
        // If e_out is above all grid points, use the last bracket
        if e_out > e_out_grid[e_out_grid.len().saturating_sub(1).min(n_ang - 1)] {
            idx = n_ang.saturating_sub(2);
        }
    }
    
    // Ensure idx is valid
    idx = idx.min(n_ang.saturating_sub(1));
    
    // Sample from angular distribution(s) at the outgoing energy
    // Interpolate between distributions if between grid points
    let mu_lo = angular_distributions[idx].sample(rng);
    let mu = if idx + 1 < n_ang && e_out_grid.len() > idx + 1 {
        let e_lo = e_out_grid[idx];
        let e_hi = e_out_grid[idx + 1];
        
        if e_hi > e_lo {
            // Sample from both distributions and interpolate
            let mu_hi = angular_distributions[idx + 1].sample(rng);
            let f = (e_out - e_lo) / (e_hi - e_lo);
            mu_lo * (1.0 - f) + mu_hi * f
        } else {
            mu_lo
        }
    } else {
        mu_lo
    };
    
    // Clamp mu to valid range
    (e_out, mu.max(-1.0).min(1.0))
}

/// Find the energy index for interpolation
fn find_energy_index(target_energy: f64, energy_grid: &[f64]) -> usize {
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



#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    
    #[test]
    fn test_correlated_single_angular_dist() {
        let energy_grid = vec![1e6];
        let energy_out = vec![TabulatedProbability::Tabulated {
            x: vec![5e5, 7e5],
            p: vec![0.5, 0.5],
        }];
        let mu_dist = vec![vec![TabulatedProbability::Tabulated {
            x: vec![-1.0, 0.0, 1.0],
            p: vec![0.33, 0.33, 0.34],
        }]];
        
        let mut rng = StdRng::seed_from_u64(42);
        let (e_out, mu) = sample_correlated_angle_energy(1e6, &energy_grid, &energy_out, &mu_dist, &mut rng);
        
        // Check valid ranges
        assert!(e_out >= 5e5 && e_out <= 7e5);
        assert!(mu >= -1.0 && mu <= 1.0);
    }
    
    #[test]
    fn test_correlated_fallback_isotropic() {
        let energy_grid = vec![1e6];
        let energy_out = vec![TabulatedProbability::Tabulated {
            x: vec![5e5],
            p: vec![1.0],
        }];
        let mu_dist = vec![vec![]]; // No angular distributions
        
        let mut rng = StdRng::seed_from_u64(42);
        let (e_out, mu) = sample_correlated_angle_energy(1e6, &energy_grid, &energy_out, &mu_dist, &mut rng);
        
        // Should fallback to isotropic
        assert!(mu >= -1.0 && mu <= 1.0);
        assert_eq!(e_out, 5e5);
    }
}
