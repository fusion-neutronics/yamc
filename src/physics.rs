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

/// Perform elastic scattering for a neutron particle
/// Uses proper two-body elastic scattering kinematics
pub fn elastic_scatter(particle: &mut Particle, awr: f64, rng: &mut impl Rng) {
    let e_in = particle.energy;

    // Sample isotropic scattering angle in center-of-mass frame
    let mu_cm = 2.0 * rng.gen::<f64>() - 1.0;

    // Calculate outgoing energy using proper elastic kinematics
    // E_out = E_in * [(AWR^2 + 1 + 2*AWR*mu_cm) / (AWR + 1)^2]
    let numerator = awr * awr + 1.0 + 2.0 * awr * mu_cm;
    let denominator = (awr + 1.0).powi(2);
    let e_out = e_in * numerator / denominator;

    particle.energy = e_out;

    // Calculate scattering angle in lab frame
    // mu_lab = (1 + AWR*mu_cm) / sqrt(1 + AWR^2 + 2*AWR*mu_cm)
    let mu_lab = (1.0 + awr * mu_cm) / (1.0 + awr * awr + 2.0 * awr * mu_cm).sqrt();

    // Sample azimuthal angle uniformly
    let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();

    // Rotate direction vector by (mu_lab, phi)
    let u_old = Vector3::from_row_slice(&particle.direction);
    let u_new = rotate_angle(u_old, mu_lab, rng);

    particle.direction = [u_new.x, u_new.y, u_new.z];
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
        };
        let awr = 1.0; // hydrogen
        elastic_scatter(&mut particle, awr, &mut rng);
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
