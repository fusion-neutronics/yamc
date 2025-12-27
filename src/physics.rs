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

/// Rotate a direction vector by angle theta (cos(theta)=mu) around arbitrary axis
/// This rotates u_old to a new direction with cosine mu relative to original
pub fn rotate_direction_3d(u_old: &Vector3<f64>, mu: f64, phi: f64) -> Vector3<f64> {
    let sin_theta = (1.0 - mu * mu).max(0.0).sqrt();

    // Find a perpendicular vector to u_old
    let perp = if u_old.x.abs() < 0.99 {
        Vector3::new(1.0, 0.0, 0.0).cross(u_old).normalize()
    } else {
        Vector3::new(0.0, 1.0, 0.0).cross(u_old).normalize()
    };
    let ortho = u_old.cross(&perp);

    mu * u_old + sin_theta * phi.cos() * perp + sin_theta * phi.sin() * ortho
}

/// Sample target velocity using the CXS (Constant Cross Section) approximation
/// This matches OpenMC's sample_cxs_target_velocity function.
/// Returns velocity in units consistent with neutron velocity (sqrt(eV))
///
/// The CXS method uses rejection sampling to properly weight the target velocity
/// distribution by the relative velocity (collision probability).
pub fn sample_cxs_target_velocity(
    awr: f64,
    neutron_energy: f64,
    neutron_direction: &[f64; 3],
    temperature_k: f64,
    rng: &mut impl Rng,
) -> Vector3<f64> {
    // Boltzmann constant in eV/K
    const K_B: f64 = 8.617333e-5;
    let k_t = K_B * temperature_k;

    // Reduced neutron velocity: beta_vn = sqrt(awr * E / kT)
    let beta_vn = (awr * neutron_energy / k_t).sqrt();

    // Probability weighting factor
    let alpha = 1.0 / (1.0 + std::f64::consts::PI.sqrt() * beta_vn / 2.0);

    let beta_vt_sq: f64;
    let mu: f64;

    loop {
        // Sample two random numbers
        let r1: f64 = rng.gen();
        let r2: f64 = rng.gen();

        let beta_vt_sq_candidate = if rng.gen::<f64>() < alpha {
            // With probability alpha, sample from p(y) = y*e^(-y)
            // Using sampling scheme: -log(r1 * r2)
            -(r1.ln() + r2.ln())
        } else {
            // With probability 1-alpha, sample from p(y) = y^2 * e^(-y^2)
            // Using sampling scheme: -log(r1) - log(r2)*cos^2(...)
            let c = (std::f64::consts::PI / 2.0 * rng.gen::<f64>()).cos();
            -r1.ln() - r2.ln() * c * c
        };

        // Determine beta * vt
        let beta_vt = beta_vt_sq_candidate.sqrt();

        // Sample cosine of angle between neutron and target velocity
        let mu_candidate = 2.0 * rng.gen::<f64>() - 1.0;

        // Determine rejection probability based on relative velocity
        // accept_prob = |v_rel| / |v_rel_max| = sqrt(vn^2 + vt^2 - 2*vn*vt*mu) / (vn + vt)
        let accept_prob = (beta_vn * beta_vn + beta_vt_sq_candidate
            - 2.0 * beta_vn * beta_vt * mu_candidate)
            .sqrt()
            / (beta_vn + beta_vt);

        // Perform rejection sampling
        if rng.gen::<f64>() < accept_prob {
            beta_vt_sq = beta_vt_sq_candidate;
            mu = mu_candidate;
            break;
        }
    }

    // Determine speed of target nucleus
    let vt = (beta_vt_sq * k_t / awr).sqrt();

    // Determine velocity vector of target nucleus based on neutron's direction
    // and the sampled angle between them
    let u = Vector3::new(
        neutron_direction[0],
        neutron_direction[1],
        neutron_direction[2],
    );

    // Rotate by angle mu around a random azimuthal angle
    let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
    vt * rotate_direction_3d(&u, mu, phi)
}

/// Legacy simple Maxwellian sampling (kept for reference/testing)
/// This does NOT properly weight by relative velocity and should not be used
/// for production simulations.
#[allow(dead_code)]
pub fn sample_target_velocity_simple(
    awr: f64,
    temperature_k: f64,
    rng: &mut impl Rng,
) -> Vector3<f64> {
    // Boltzmann constant in eV/K
    const K_B: f64 = 8.617333e-5;

    // Thermal energy of target nucleus
    let k_t = K_B * temperature_k;

    // Mass of target in neutron masses
    let mass_target = awr;

    // Sample velocity components from Gaussian distribution
    // v ~ N(0, sqrt(kT/m))
    use rand_distr::{Distribution, Normal};
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

    // Sample target velocity using CXS approximation (matches OpenMC)
    let v_t = sample_cxs_target_velocity(awr, particle.energy, &particle.direction, temperature_k, rng);

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
//   FISSION PHYSICS
// =====================

/// Sample fission neutron energy from Watt spectrum
/// Uses the standard Watt fission spectrum: χ(E) = C * exp(-E/a) * sinh(sqrt(b*E))
/// Default parameters are for U-235 thermal fission: a = 0.988 MeV, b = 2.249 MeV^-1
/// Returns energy in eV
pub fn sample_watt_spectrum(rng: &mut impl Rng) -> f64 {
    // Watt spectrum parameters for U-235 (default, reasonable for actinides)
    // a and b in MeV
    let a = 0.988;
    let b = 2.249;

    // Rejection sampling for Watt spectrum
    // Using the method from LA-UR-14-27694 (OpenMC theory manual)
    let k = 1.0 + (b * a / 8.0);

    loop {
        let r1: f64 = rng.gen();
        let r2: f64 = rng.gen();
        let r3: f64 = rng.gen();

        // Sample from g(E) = C * exp(-(E/a - sqrt(b*a/4))^2)
        let w = a * k + (-(a * k * r1).ln()) * (r2.cos().powi(2));

        // Rejection test
        let eta = (b * w).sqrt();
        if r3 <= (eta.sinh() / eta) {
            // Convert from MeV to eV
            return w * 1.0e6;
        }
    }
}

/// Sample fission neutrons and return them as a vector of particles
///
/// # Arguments
/// * `particle` - The incident neutron particle
/// * `nu_bar` - Average number of neutrons per fission at this energy
/// * `fission_products` - Optional fission neutron products with energy distributions
/// * `rng` - Random number generator
///
/// # Returns
/// Vector of fission neutron particles (may be empty in rare cases)
pub fn sample_fission_neutrons<R: Rng>(
    particle: &crate::particle::Particle,
    nu_bar: f64,
    fission_products: Option<&[&crate::reaction_product::ReactionProduct]>,
    rng: &mut R,
) -> Vec<crate::particle::Particle> {
    // Stochastic rounding to determine actual number of neutrons
    let n_neutrons = {
        let base = nu_bar.floor() as usize;
        let frac = nu_bar - nu_bar.floor();
        if rng.gen::<f64>() < frac {
            base + 1
        } else {
            base
        }
    };

    if n_neutrons == 0 {
        return Vec::new();
    }

    let mut neutrons = Vec::with_capacity(n_neutrons);

    // Check if we have product distributions to sample from
    let use_products = fission_products
        .map(|prods| !prods.is_empty())
        .unwrap_or(false);

    for _ in 0..n_neutrons {
        let mut new_particle = particle.clone();

        if use_products {
            // Sample from product distributions (usually prompt neutron spectrum)
            let products = fission_products.unwrap();
            let prompt_product = products[0]; // First product is typically prompt neutrons

            let (e_fission, mu) = prompt_product.sample(particle.energy, rng);
            new_particle.energy = e_fission.max(1e-11); // Ensure positive energy

            // Rotate direction by sampled mu
            crate::inelastic::rotate_direction(&mut new_particle.direction, mu, rng);
        } else {
            // No product data - use Watt spectrum with isotropic emission
            new_particle.energy = sample_watt_spectrum(rng);

            // Isotropic direction
            let mu = 2.0 * rng.gen::<f64>() - 1.0;
            let phi = 2.0 * std::f64::consts::PI * rng.gen::<f64>();
            let sin_theta = (1.0 - mu * mu).sqrt();
            new_particle.direction = [
                sin_theta * phi.cos(),
                sin_theta * phi.sin(),
                mu,
            ];
        }

        neutrons.push(new_particle);
    }

    neutrons
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
