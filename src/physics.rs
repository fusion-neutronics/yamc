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

/// Perform elastic scattering for a neutron particle with thermal motion
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
