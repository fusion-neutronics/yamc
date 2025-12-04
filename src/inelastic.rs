use crate::particle::Particle;
use crate::reaction::Reaction;
use crate::reaction_product::ParticleType;
use rand::Rng;

/// Determine neutron multiplicity for inelastic scattering based on MT number
/// Returns the number of neutrons produced for inelastic reactions only
/// Does not handle fission or other non-inelastic reactions
pub fn get_inelastic_neutron_multiplicity(mt: i32) -> usize {
    match mt {
        // Inelastic scattering levels (MT 51-90) - 1 neutron out
        51..=90 => 1,
        
        // Continuum inelastic - 1 neutron out
        91 => 1,
        
        // (n,2n) reactions - 2 neutrons out
        16 => 2,
        
        // (n,3n) reactions - 3 neutrons out  
        17 => 3,
        
        // (n,4n) reactions - 4 neutrons out
        37 => 4,
        
        // (n,n') first excited state - 1 neutron out
        4 => 1,
        
        // Other inelastic reactions that produce neutrons
        // MT 22: (n,n'alpha) - 1 neutron + alpha
        22 => 1,
        
        // MT 28: (n,n'p) - 1 neutron + proton  
        28 => 1,
        
        // MT 32: (n,n'd) - 1 neutron + deuteron
        32 => 1,
        
        // MT 33: (n,n't) - 1 neutron + triton
        33 => 1,
        
        // MT 34: (n,n'3He) - 1 neutron + 3He
        34 => 1,
        
        // Panic for unknown MT numbers to catch unsupported reaction types
        _ => panic!("Unsupported inelastic reaction MT number: {}", mt),
    }
}

/// Handle inelastic scattering following OpenMC's approach
/// Returns vector of outgoing neutrons based on reaction products or analytical models
pub fn inelastic_scatter<R: rand::Rng>(
    particle: &Particle,
    reaction: &Reaction,
    nuclide_name: &str, // nuclide name for AWR lookup when needed
    rng: &mut R,
) -> Vec<Particle> {
    // Check if reaction has product data
    if !reaction.products.is_empty() {
        // Sample from product distributions when available
        sample_from_products(particle, reaction, rng)
    } else {
        // No product data available - use analytical approach based on Q-value and MT
        analytical_inelastic_scatter(particle, reaction, nuclide_name, rng)
    }
}

/// Sample outgoing particles from reaction product distributions
/// This handles the case where explicit product data is available
/// This is public so it can be used for MT=2 elastic scattering with product data
pub fn sample_from_products<R: rand::Rng>(
    particle: &Particle,
    reaction: &Reaction,
    rng: &mut R,
) -> Vec<Particle> {
    let incoming_energy = particle.energy;
    
    // Find neutron products (since we're tracking neutrons)
    let neutron_products: Vec<&crate::reaction_product::ReactionProduct> = reaction
        .products
        .iter()
        .filter(|product| product.is_particle_type(&ParticleType::Neutron))
        .collect();
    
    if neutron_products.is_empty() {
        // No neutron products - particle is absorbed (e.g., (n,gamma), (n,p), etc.)
        return Vec::new();
    }
    
    // Create outgoing neutrons for each neutron product
    let mut outgoing_neutrons = Vec::new();
    
    for neutron_product in neutron_products {
        // Sample outgoing energy and scattering cosine from product distribution
        let (e_out, mu) = neutron_product.sample(incoming_energy, rng);
        
        // Create new neutron particle
        let mut new_particle = particle.clone();
        new_particle.energy = e_out;
        rotate_direction(&mut new_particle.direction, mu, rng);
        
        outgoing_neutrons.push(new_particle);
    }
    
    outgoing_neutrons
}

/// Analytical inelastic scattering based on Q-values and reaction kinematics
/// Follows OpenMC's LevelInelastic approach for MT 51-90 and general inelastic for others
/// This handles reactions without product data using analytical energy-angle relationships
pub fn analytical_inelastic_scatter<R: rand::Rng>(
    particle: &Particle,
    reaction: &Reaction,
    nuclide_name: &str,
    rng: &mut R,
) -> Vec<Particle> {
    use crate::data::ATOMIC_WEIGHT_RATIO;
    
    let e_in = particle.energy;
    let mt = reaction.mt_number;
    let q_value = reaction.q_value;
    
    // Determine neutron multiplicity based on reaction type
    // TODO: Replace with actual yield data from nuclear data files when available
    let neutron_multiplicity = get_inelastic_neutron_multiplicity(mt);
    
    // Look up atomic weight ratio only when needed for analytical calculations
    let awr = *ATOMIC_WEIGHT_RATIO
        .get(nuclide_name)
        .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));

    let e_out = if mt >= 51 && mt <= 90 {
        // Level inelastic scattering (MT 51-90) - OpenMC LevelInelastic approach
        // E_out = mass_ratio * (E_in - threshold) where threshold = (A+1)/A * |Q|
        let threshold = (awr + 1.0) / awr * q_value.abs();
        let mass_ratio = (awr / (awr + 1.0)).powi(2);
        
        if e_in < threshold {
            eprintln!("Warning: Energy {} eV below threshold {} eV for MT {}", e_in, threshold, mt);
            return Vec::new();
        }
        
        mass_ratio * (e_in - threshold)
    } else if reaction.products.is_empty() {
        // No product data available - use Q-value based analytical approach (OpenMC fallback)
        if q_value < 0.0 {
            // Endothermic reaction
            let threshold = (awr + 1.0) / awr * q_value.abs();
            if e_in < threshold {
                eprintln!("Warning: Energy {} eV below endothermic threshold {} eV for MT {}", e_in, threshold, mt);
                return Vec::new();
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

    // For multi-neutron reactions, divide energy among neutrons
    let energy_per_neutron = if neutron_multiplicity > 1 {
        // Simple energy sharing model - could be more sophisticated
        e_out / neutron_multiplicity as f64
    } else {
        e_out
    };

    // Create the specified number of outgoing neutrons
    let mut outgoing_neutrons = Vec::with_capacity(neutron_multiplicity);
    
    for _ in 0..neutron_multiplicity {
        // Create new neutron particle
        let mut new_particle = particle.clone();
        new_particle.energy = energy_per_neutron;
        
        // Isotropic scattering angle (OpenMC default for missing angular data)
        let mu = rng.gen_range(-1.0..=1.0);
        rotate_direction(&mut new_particle.direction, mu, rng);
        
        outgoing_neutrons.push(new_particle);
    }
    
    outgoing_neutrons
}

/// Rotate particle direction by scattering angle with cosine mu
/// This implements isotropic azimuthal angle sampling
pub fn rotate_direction<R: rand::Rng>(
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