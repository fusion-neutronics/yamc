use crate::particle::Particle;
use crate::reaction::Reaction;
use crate::reaction_product::ParticleType;
use rand::Rng;

/// Handle inelastic scattering following OpenMC's approach
/// Updates particle energy and direction based on reaction products or analytical models
pub fn inelastic_scatter<R: rand::Rng>(
    particle: &mut Particle,
    reaction: &Reaction,
    awr: f64, // atomic weight ratio
    rng: &mut R,
) {
    // Check if reaction has product data
    if !reaction.products.is_empty() {
        // Sample from product distributions when available
        sample_from_products(particle, reaction, rng);
    } else {
        // No product data available - use analytical approach based on Q-value and MT
        analytical_inelastic_scatter(particle, reaction, awr, rng);
    }
}

/// Sample outgoing particle from reaction product distributions
/// This handles the case where explicit product data is available
fn sample_from_products<R: rand::Rng>(
    particle: &mut Particle,
    reaction: &Reaction,
    rng: &mut R,
) {
    let incoming_energy = particle.energy;
    
    // Find neutron products (since we're tracking neutrons)
    let neutron_products: Vec<&crate::reaction_product::ReactionProduct> = reaction
        .products
        .iter()
        .filter(|product| product.is_particle_type(&ParticleType::Neutron))
        .collect();
    
    if neutron_products.is_empty() {
        panic!("Reaction has products but no neutron products found - data inconsistency for MT {}", reaction.mt_number);
    }
    
    // For inelastic scattering, we typically expect one neutron product
    // If multiple neutron products, sample from the first one (most common case)
    // More sophisticated sampling could handle multiple products
    let neutron_product = neutron_products[0];
    
    // Sample outgoing energy and scattering cosine from product distribution
    let (e_out, mu) = neutron_product.sample(incoming_energy, rng);
    
    // Update particle energy and direction
    particle.energy = e_out;
    rotate_direction(&mut particle.direction, mu, rng);
}

/// Analytical inelastic scattering based on Q-values and reaction kinematics
/// Follows OpenMC's LevelInelastic approach for MT 51-90 and general inelastic for others
/// This handles reactions without product data using analytical energy-angle relationships
fn analytical_inelastic_scatter<R: rand::Rng>(
    particle: &mut Particle,
    reaction: &Reaction,
    awr: f64,
    rng: &mut R,
) {
    let e_in = particle.energy;
    let mt = reaction.mt_number;
    let q_value = reaction.q_value;

    let e_out = if mt >= 51 && mt <= 90 {
        // Level inelastic scattering (MT 51-90) - OpenMC LevelInelastic approach
        // E_out = mass_ratio * (E_in - threshold) where threshold = (A+1)/A * |Q|
        let threshold = (awr + 1.0) / awr * q_value.abs();
        let mass_ratio = (awr / (awr + 1.0)).powi(2);
        
        if e_in < threshold {
            eprintln!("Warning: Energy {} eV below threshold {} eV for MT {}", e_in, threshold, mt);
            return;
        }
        
        mass_ratio * (e_in - threshold)
    } else if reaction.products.is_empty() {
        // No product data available - use Q-value based analytical approach (OpenMC fallback)
        if q_value < 0.0 {
            // Endothermic reaction
            let threshold = (awr + 1.0) / awr * q_value.abs();
            if e_in < threshold {
                eprintln!("Warning: Energy {} eV below endothermic threshold {} eV for MT {}", e_in, threshold, mt);
                return;
            }
            e_in + q_value * awr / (awr + 1.0)
        } else {
            // Exothermic reaction
            e_in + q_value * awr / (awr + 1.0)
        }
    } else {
        // Has products but using analytical fallback - simplified energy model
        let max_energy_loss = if q_value < 0.0 {
            (-q_value * (awr + 1.0) / awr).min(e_in * 0.9)
        } else {
            0.0
        };
        
        if max_energy_loss > 0.0 {
            let energy_loss = rng.gen_range(0.0..max_energy_loss);
            (e_in - energy_loss).max(0.1)
        } else {
            e_in + q_value.max(0.0)
        }
    };

    // Isotropic scattering angle (OpenMC default for missing angular data)
    let mu = rng.gen_range(-1.0..=1.0);
    particle.energy = e_out;
    rotate_direction(&mut particle.direction, mu, rng);
}

/// Rotate particle direction by scattering angle with cosine mu
/// This implements isotropic azimuthal angle sampling
fn rotate_direction<R: rand::Rng>(
    direction: &mut [f64; 3],
    mu: f64, // cosine of scattering angle
    rng: &mut R,
) {
    // Sample azimuthal angle uniformly
    let phi = rng.gen_range(0.0..2.0 * std::f64::consts::PI);
    
    // Current direction
    let [ux, uy, uz] = *direction;
    
    // Calculate sine of scattering angle
    let sin_theta = (1.0 - mu * mu).sqrt();
    
    // For isotropic scattering, we need to handle the case where uz ≈ ±1
    let (new_ux, new_uy, new_uz) = if uz.abs() < 0.999 {
        // General case
        let factor = sin_theta / (1.0 - uz * uz).sqrt();
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        
        (
            mu * ux + factor * (ux * uz * cos_phi - uy * sin_phi),
            mu * uy + factor * (uy * uz * cos_phi + ux * sin_phi),
            mu * uz - factor * (1.0 - uz * uz).sqrt() * cos_phi,
        )
    } else {
        // Special case: uz ≈ ±1 (beam nearly parallel to z-axis)
        let cos_phi = phi.cos();
        let sin_phi = phi.sin();
        let sign = if uz > 0.0 { 1.0 } else { -1.0 };
        
        (
            sin_theta * cos_phi,
            sin_theta * sin_phi,
            sign * mu,
        )
    };
    
    *direction = [new_ux, new_uy, new_uz];
    
    // Normalize to ensure unit vector (numerical precision)
    let norm = (new_ux * new_ux + new_uy * new_uy + new_uz * new_uz).sqrt();
    if norm > 1e-10 {
        direction[0] /= norm;
        direction[1] /= norm;
        direction[2] /= norm;
    }
}