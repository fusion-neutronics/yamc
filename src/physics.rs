// Elastic scattering physics for Monte Carlo transport

use crate::particle::Particle;
use nalgebra::Vector3;
use rand::Rng;

/// Rotate a vector by a random angle mu_cm around a random axis
pub fn rotate_angle(u_cm: Vector3<f64>, mu_cm: f64, rng: &mut impl Rng) -> Vector3<f64> {
    let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
    let perp = if u_cm.x.abs() < 0.99 {
        Vector3::new(1.0, 0.0, 0.0).cross(&u_cm).normalize()
    } else {
        Vector3::new(0.0, 1.0, 0.0).cross(&u_cm).normalize()
    };
    let ortho = u_cm.cross(&perp);
    let sin_theta = (1.0 - mu_cm * mu_cm).sqrt();
    mu_cm * u_cm + sin_theta * phi.cos() * perp + sin_theta * phi.sin() * ortho
}

/// Sample target velocity from Maxwell-Boltzmann distribution
/// Returns velocity in units consistent with neutron velocity (sqrt(eV))
fn sample_target_velocity(awr: f64, temperature_k: f64, rng: &mut impl Rng) -> Vector3<f64> {
    // Boltzmann constant in eV/K
    const K_B: f64 = 8.617333e-5;
    
    // Thermal energy of target nucleus
    let k_t = K_B * temperature_k;
    
    // Mass of target in neutron masses
    let mass_target = awr;
    
    // Sample velocity components from Gaussian distribution
    // v ~ N(0, sqrt(kT/m))
    use rand_distr::{Normal, Distribution};
    let sigma = (k_t / mass_target).sqrt();
    let normal = Normal::new(0.0, sigma).unwrap();
    
    Vector3::new(
        normal.sample(rng),
        normal.sample(rng),
        normal.sample(rng),
    )
}

/// Perform elastic scattering for a neutron particle with thermal motion (isotropic in CM)
/// This is used when no angular distribution data is available
pub fn elastic_scatter(particle: &mut Particle, awr: f64, temperature_k: f64, rng: &mut impl Rng) {
    // Neutron velocity in LAB (velocity = sqrt(energy) in our units where neutron mass = 1)
    let vel = particle.energy.sqrt();
    let v_n = Vector3::from_row_slice(&particle.direction) * vel;

    // Sample target velocity from Maxwell-Boltzmann distribution
    let v_t = sample_target_velocity(awr, temperature_k, rng);

    // Center-of-mass velocity
    let v_cm = (v_n + awr * v_t) / (awr + 1.0);

    // Transform to CM frame
    let mut v_n_cm = v_n - v_cm;
    let vel_cm = v_n_cm.norm();
    
    // Handle case where neutron and target have same velocity
    if vel_cm < 1e-10 {
        // Just scatter isotropically in lab frame
        let mu = 2.0 * rng.gen::<f64>() - 1.0;
        let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
        let sin_theta = (1.0 - mu * mu).sqrt();
        particle.direction = [
            sin_theta * phi.cos(),
            sin_theta * phi.sin(),
            mu,
        ];
        return;
    }

    // Sample scattering angle in CM (isotropic for elastic)
    let mu_cm = 2.0 * rng.gen::<f64>() - 1.0;

    // Direction in CM
    let u_cm = v_n_cm / vel_cm;

    // Rotate neutron velocity vector to new angle
    v_n_cm = vel_cm * rotate_angle(u_cm, mu_cm, rng);

    // Transform back to LAB frame
    let v_n_lab = v_n_cm + v_cm;

    // Update particle energy and direction
    particle.energy = v_n_lab.dot(&v_n_lab);
    let vel_lab = particle.energy.sqrt();
    
    if vel_lab > 1e-10 {
        particle.direction = [
            v_n_lab.x / vel_lab,
            v_n_lab.y / vel_lab,
            v_n_lab.z / vel_lab,
        ];
    }
}

/// DBRC (Doppler Broadening Rejection Correction) elastic scattering
/// Uses target-at-rest approximation with tabulated angular distribution
/// This is appropriate for the resonance region (typically ~1 eV to 210 keV)
pub fn dbrc_elastic_scatter(
    particle: &mut Particle, 
    awr: f64, 
    // temperature_k: f64,
    mu_cm: f64,  // Scattering cosine from tabulated distribution
    rng: &mut impl Rng
) {
    // For DBRC, OpenMC uses target-at-rest kinematics with the tabulated angular distribution
    // The thermal motion is already accounted for in the cross section sampling via relative energy
    
    // Target-at-rest: transform using simple kinematics
    // E_out = E_in * ((A^2 + 1 + 2*A*mu) / (A+1)^2) for elastic scattering
    
    let e_in = particle.energy;
    let a = awr;
    
    // Scattering cosine in lab frame (from CM using target-at-rest kinematics)
    // mu_lab = (1 + A*mu_cm) / sqrt(A^2 + 2*A*mu_cm + 1)
    let mu_lab = (1.0 + a * mu_cm) / (a * a + 2.0 * a * mu_cm + 1.0).sqrt();
    
    // Outgoing energy (target-at-rest kinematics)
    let e_out = e_in * (a * a + 1.0 + 2.0 * a * mu_cm) / ((a + 1.0) * (a + 1.0));
    
    // Update particle energy
    particle.energy = e_out;
    
    // Rotate direction
    let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
    let sin_theta = (1.0 - mu_lab * mu_lab).max(0.0).sqrt();
    
    // Create perpendicular vectors to current direction
    let old_dir = Vector3::from_row_slice(&particle.direction);
    let perp = if old_dir.x.abs() < 0.99 {
        Vector3::new(1.0, 0.0, 0.0).cross(&old_dir).normalize()
    } else {
        Vector3::new(0.0, 1.0, 0.0).cross(&old_dir).normalize()
    };
    let ortho = old_dir.cross(&perp);
    
    // New direction
    let new_dir = mu_lab * old_dir + sin_theta * phi.cos() * perp + sin_theta * phi.sin() * ortho;
    
    particle.direction = [new_dir.x, new_dir.y, new_dir.z];
}

// =====================
//        TESTS
// =====================

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::SeedableRng;

    #[test]
    fn test_rotate_angle_preserves_norm() {
        let mut rng = StdRng::seed_from_u64(42);
        let u_cm = Vector3::new(0.0, 0.0, 1.0);
        let mu_cm = 0.5;
        let v = rotate_angle(u_cm, mu_cm, &mut rng);
        // Should be unit vector
        assert!((v.norm() - 1.0).abs() < 1e-12, "norm = {}", v.norm());
        // z-component should be mu_cm
        assert!((v.z - mu_cm).abs() < 1e-12, "z = {} mu_cm = {}", v.z, mu_cm);
    }

    #[test]
    fn test_elastic_scatter_energy_and_direction() {
        let mut rng = StdRng::seed_from_u64(123);
        // Dummy particle: energy = 2.0, direction = [0,0,1]
        let mut particle = Particle {
            energy: 2.0,
            direction: [0.0, 0.0, 1.0],
            alive: true,
            position: [0.0, 0.0, 0.0],
            id: 0,
            current_cell_index: None,
        };
        let awr = 1.0; // hydrogen
        let temperature_k = 294.0; // room temperature
        elastic_scatter(&mut particle, awr, temperature_k, &mut rng);
        // Energy should be positive
        assert!(particle.energy > 0.0);
        // Direction should be unit vector
        let norm = (particle.direction[0].powi(2)
            + particle.direction[1].powi(2)
            + particle.direction[2].powi(2))
        .sqrt();
        assert!((norm - 1.0).abs() < 1e-12, "norm = {}", norm);
    }
}
