// Kalbach-Mann correlated angle-energy distribution
// Corresponds to OpenMC's secondary_kalbach.cpp

use rand::Rng;
use crate::reaction_product::{Tabulated1D, TabulatedProbability};

/// Sample from Kalbach-Mann correlated angle-energy distribution
/// 
/// This function corresponds to OpenMC's KalbachMann::sample()
/// 
/// The Kalbach-Mann systematics provides a correlation between outgoing energy
/// and scattering angle for pre-equilibrium and direct reactions. The angular
/// distribution is given by: f(mu) ∝ exp(a*mu) where 'a' is the Kalbach slope
/// parameter that depends on outgoing energy.
/// 
/// # Arguments
/// * `incoming_energy` - Incident particle energy (eV)
/// * `energy_grid` - Tabulated incident energy points
/// * `energy_out` - Outgoing energy probability distributions at each incident energy
/// * `slope` - Kalbach slope parameter 'a' as function of outgoing energy
/// * `rng` - Random number generator
/// 
/// # Returns
/// Tuple of (outgoing_energy, mu_cosine)
pub fn sample_kalbach_mann<R: Rng>(
    incoming_energy: f64,
    energy_grid: &[f64],
    energy_out: &[TabulatedProbability],
    slope: &[Tabulated1D],
    rng: &mut R,
) -> (f64, f64) {
    // Find energy bracket for incoming energy
    let i = find_energy_index(incoming_energy, energy_grid);
    
    if i >= energy_out.len() || i >= slope.len() {
        // Fallback to isotropic if data missing
        return (incoming_energy, 2.0 * rng.gen::<f64>() - 1.0);
    }
    
    // Sample outgoing energy from tabulated distribution
    let e_out = energy_out[i].sample(rng);
    
    // Evaluate Kalbach slope parameter at sampled outgoing energy
    let a = slope[i].evaluate(e_out);
    
    // Sample scattering angle from Kalbach-Mann angular distribution
    // f(mu) = exp(a*mu) / [2*sinh(a)/a] for a != 0
    // f(mu) = 1/2 for a = 0 (isotropic)
    let mu = if a.abs() < 1e-6 {
        // Nearly isotropic case
        2.0 * rng.gen::<f64>() - 1.0
    } else {
        // Sample from exp(a*mu) on [-1, 1]
        // Using inverse CDF method:
        // CDF(mu) = [exp(a*mu) - exp(-a)] / [exp(a) - exp(-a)]
        // Solving for mu: mu = ln[(1-r)*exp(-a) + r*exp(a)] / a
        let r = rng.gen::<f64>();
        if a > 0.0 {
            ((1.0 - r) * (-a).exp() + r * a.exp()).ln() / a
        } else {
            // For negative a, rearrange to avoid numerical issues
            -((1.0 - r) * a.exp() + r * (-a).exp()).ln() / a
        }
    };
    
    // Clamp mu to valid range [-1, 1]
    let mu = mu.max(-1.0).min(1.0);
    
    (e_out, mu)
}

/// Find the energy index for interpolation
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

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;
    
    #[test]
    fn test_kalbach_mann_isotropic() {
        // Test isotropic case (a ≈ 0)
        let energy_grid = vec![1e6];
        let energy_out = vec![TabulatedProbability::Tabulated {
            x: vec![5e5, 7e5],
            p: vec![0.5, 0.5],
        }];
        let slope = vec![crate::reaction_product::Tabulated1D::Tabulated1D {
            x: vec![5e5, 7e5],
            y: vec![0.0, 0.0],  // a = 0 => isotropic
            breakpoints: vec![],
            interpolation: vec![],
        }];
        
        let mut rng = StdRng::seed_from_u64(42);
        let (e_out, mu) = sample_kalbach_mann(1e6, &energy_grid, &energy_out, &slope, &mut rng);
        
        // Check mu is in valid range
        assert!(mu >= -1.0 && mu <= 1.0);
        // Check energy was sampled
        assert!(e_out >= 5e5 && e_out <= 7e5);
    }
    
    #[test]
    fn test_kalbach_mann_forward_peaked() {
        // Test forward-peaked case (a > 0)
        let energy_grid = vec![1e6];
        let energy_out = vec![TabulatedProbability::Tabulated {
            x: vec![5e5],
            p: vec![1.0],
        }];
        let slope = vec![crate::reaction_product::Tabulated1D::Tabulated1D {
            x: vec![5e5],
            y: vec![5.0],  // a = 5.0 => strongly forward peaked
            breakpoints: vec![],
            interpolation: vec![],
        }];
        
        let mut rng = StdRng::seed_from_u64(42);
        
        // Sample many times and check that mu is mostly positive (forward)
        let mut forward_count = 0;
        for _ in 0..100 {
            let (_, mu) = sample_kalbach_mann(1e6, &energy_grid, &energy_out, &slope, &mut rng);
            assert!(mu >= -1.0 && mu <= 1.0);
            if mu > 0.0 {
                forward_count += 1;
            }
        }
        
        // With a = 5.0, most samples should be forward (mu > 0)
        assert!(forward_count > 70, "Forward peaked distribution should favor mu > 0");
    }
}
