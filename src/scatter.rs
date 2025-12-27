use crate::particle::Particle;
use crate::reaction::Reaction;
use crate::reaction_product::ParticleType;

/// Determine neutron multiplicity for scattering reactions based on MT number
/// These are scattering reactions that produce neutrons (excluding MT 2 elastic and MT 50-91 inelastic)
/// Returns the number of neutrons produced
pub fn get_scatter_neutron_multiplicity(mt: i32) -> usize {
    match mt {
        // Single neutron producing reactions
        5 => 1,   // (n,anything) - usually 1 neutron
        11 => 1,  // (n,2nd) - 1 neutron + deuteron
        22 => 1,  // (n,n'alpha) - 1 neutron + alpha
        23 => 1,  // (n,n'3alpha) - 1 neutron + 3 alpha
        24 => 1,  // (n,2n'alpha) - treated as 1n for simplicity
        25 => 1,  // (n,3n'alpha) - treated as 1n for simplicity
        28 => 1,  // (n,n'p) - 1 neutron + proton
        29 => 1,  // (n,n'2alpha) - 1 neutron + 2 alpha
        30 => 1,  // (n,2n'2alpha) - treated as 1n for simplicity
        32 => 1,  // (n,n'd) - 1 neutron + deuteron
        33 => 1,  // (n,n't) - 1 neutron + triton
        34 => 1,  // (n,n'He-3) - 1 neutron + He-3
        35 => 1,  // (n,n'd2alpha) - 1 neutron + deuteron + 2 alpha
        36 => 1,  // (n,n't2alpha) - 1 neutron + triton + 2 alpha
        41 => 1,  // (n,2n'p) - treated as 1n for simplicity
        42 => 1,  // (n,3n'p) - treated as 1n for simplicity
        44 => 1,  // (n,n'2p) - 1 neutron + 2 protons
        45 => 1,  // (n,n'p'alpha) - 1 neutron + proton + alpha

        // Multi-neutron producing reactions (actually output multiple neutrons)
        16 => 2,  // (n,2n) - 2 neutrons out
        17 => 3,  // (n,3n) - 3 neutrons out
        37 => 4,  // (n,4n) - 4 neutrons out

        // MT 152-200: Various complex neutron-producing reactions
        152 => 5, // (n,5n) - 5 neutrons out
        153 => 1, // Complex reaction - 1 neutron
        154 => 1, // (n,2n't) or similar - 1 neutron
        156 => 1, // (n,2n'He-3) or similar - 1 neutron
        157 => 1, // (n,3n'd) or similar - 1 neutron
        158 => 1, // (n,3n't) or similar - 1 neutron
        159 => 1, // (n,3n'He-3) or similar - 1 neutron
        160 => 1, // (n,4n'p) or similar - 1 neutron
        161 => 1, // (n,4n'd) or similar - 1 neutron
        162 => 1, // (n,4n't) or similar - 1 neutron
        163 => 1, // (n,4n'He-3) or similar - 1 neutron
        164 => 1, // (n,4n'alpha) or similar - 1 neutron
        165 => 1, // (n,5n'p) or similar - 1 neutron
        166 => 1, // (n,5n'd) or similar - 1 neutron
        167 => 1, // (n,5n't) or similar - 1 neutron
        168 => 1, // (n,5n'He-3) or similar - 1 neutron
        169 => 1, // (n,5n'alpha) or similar - 1 neutron
        170 => 1, // (n,6n'p) or similar - 1 neutron
        171 => 1, // (n,6n'd) or similar - 1 neutron
        172 => 1, // (n,6n't) or similar - 1 neutron
        173 => 1, // (n,6n'He-3) or similar - 1 neutron
        174 => 1, // (n,6n'alpha) or similar - 1 neutron
        175 => 1, // (n,7n'p) or similar - 1 neutron
        176 => 1, // (n,7n'd) or similar - 1 neutron
        177 => 1, // (n,7n't) or similar - 1 neutron
        178 => 1, // (n,7n'He-3) or similar - 1 neutron
        179 => 1, // (n,7n'alpha) or similar - 1 neutron
        180 => 1, // (n,8n'p) or similar - 1 neutron
        181 => 1, // (n,8n'd) or similar - 1 neutron
        183 => 1, // (n,8n'He-3) or similar - 1 neutron
        184 => 1, // (n,8n'alpha) or similar - 1 neutron
        185 => 1, // (n,9n'p) or similar - 1 neutron
        186 => 1, // (n,9n'd) or similar - 1 neutron
        187 => 1, // (n,9n't) or similar - 1 neutron
        188 => 1, // (n,9n'He-3) or similar - 1 neutron
        189 => 1, // (n,9n'alpha) or similar - 1 neutron
        190 => 1, // (n,10n'p) or similar - 1 neutron
        194 => 5, // (n,5n) - 5 neutrons out
        195 => 6, // (n,6n) - 6 neutrons out
        196 => 7, // (n,7n) - 7 neutrons out
        198 => 8, // (n,8n) - 8 neutrons out
        199 => 9, // (n,9n) - 9 neutrons out
        200 => 10, // (n,10n) - 10 neutrons out

        // Panic for unknown MT numbers
        _ => panic!("Unsupported scattering reaction MT number: {}", mt),
    }
}

/// Handle general scattering reactions (non-elastic, non-inelastic-level)
/// This includes reactions like (n,2n), (n,3n), (n,n'alpha), etc.
/// Returns vector of outgoing neutrons based on reaction products or analytical models
pub fn scatter<R: rand::Rng>(
    particle: &Particle,
    reaction: &Reaction,
    nuclide_name: &str, // nuclide name for AWR lookup when needed
    rng: &mut R,
) -> Vec<Particle> {
    use crate::data::ATOMIC_WEIGHT_RATIO;
    let awr = *ATOMIC_WEIGHT_RATIO
        .get(nuclide_name)
        .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));
    scatter_with_awr(particle, reaction, nuclide_name, awr, rng)
}

/// Handle general scattering reactions with explicit AWR for CM to LAB conversion
pub fn scatter_with_awr<R: rand::Rng>(
    particle: &Particle,
    reaction: &Reaction,
    nuclide_name: &str,
    awr: f64,
    rng: &mut R,
) -> Vec<Particle> {
    // Check if reaction has product data
    if !reaction.products.is_empty() {
        // Sample from product distributions when available (with AWR for CM to LAB)
        sample_from_products_with_awr(particle, reaction, awr, rng)
    } else {
        // No product data available - use analytical approach based on Q-value and MT
        analytical_scatter(particle, reaction, nuclide_name, rng)
    }
}

/// Analytical scattering based on Q-values and reaction kinematics
/// Similar to inelastic scatter but uses scatter-specific multiplicity
fn analytical_scatter<R: rand::Rng>(
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
    let neutron_multiplicity = get_scatter_neutron_multiplicity(mt);

    // Look up atomic weight ratio
    let awr = *ATOMIC_WEIGHT_RATIO
        .get(nuclide_name)
        .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));

    // Calculate outgoing energy using Q-value
    let e_out = if q_value < 0.0 {
        // Endothermic reaction
        let threshold = (awr + 1.0) / awr * q_value.abs();
        if e_in < threshold {
            eprintln!("Warning: Energy {} eV below threshold {} eV for MT {}", e_in, threshold, mt);
            return Vec::new();
        }
        e_in + q_value * awr / (awr + 1.0)
    } else {
        // Exothermic reaction
        e_in + q_value * awr / (awr + 1.0)
    };

    // For multi-neutron reactions, divide energy among neutrons
    let energy_per_neutron = if neutron_multiplicity > 1 {
        // Simple energy sharing model
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

        // Isotropic scattering angle (default for missing angular data)
        let mu = rng.gen_range(-1.0..=1.0);
        crate::inelastic::rotate_direction(&mut new_particle.direction, mu, rng);

        outgoing_neutrons.push(new_particle);
    }

    outgoing_neutrons
}

/// Sample outgoing particles from reaction product distributions
/// This handles the case where explicit product data is available
fn sample_from_products<R: rand::Rng>(
    particle: &Particle,
    reaction: &Reaction,
    rng: &mut R,
) -> Vec<Particle> {
    sample_from_products_with_awr(particle, reaction, 1.0, rng)
}

/// Sample outgoing particles with explicit AWR for CM to LAB conversion
/// Follows OpenMC's approach: evaluate yield and create appropriate number of particles
pub fn sample_from_products_with_awr<R: rand::Rng>(
    particle: &Particle,
    reaction: &Reaction,
    awr: f64,
    rng: &mut R,
) -> Vec<Particle> {
    let e_in = particle.energy;

    // Find neutron products (since we're tracking neutrons)
    let neutron_products: Vec<&crate::reaction_product::ReactionProduct> = reaction
        .products
        .iter()
        .filter(|product| product.is_particle_type(&ParticleType::Neutron))
        .collect();

    if neutron_products.is_empty() {
        // No neutron products - particle is absorbed (shouldn't happen for scattering)
        return Vec::new();
    }

    // Create outgoing neutrons for each neutron product
    let mut outgoing_neutrons = Vec::new();

    for neutron_product in neutron_products {
        // Evaluate yield at incident energy (OpenMC: physics.cpp line 1157)
        let yield_val = neutron_product
            .product_yield
            .as_ref()
            .map(|y| y.evaluate(e_in))
            .unwrap_or(1.0);

        // Determine number of particles to create based on yield
        // OpenMC: if yield is integral, create that many particles
        // If non-integral, use stochastic rounding
        let n_particles = if (yield_val - yield_val.floor()).abs() < 1e-10 {
            // Integral yield - create exactly that many
            yield_val.round() as usize
        } else {
            // Non-integral yield - stochastic rounding
            let base = yield_val.floor() as usize;
            let frac = yield_val - yield_val.floor();
            if rng.gen::<f64>() < frac {
                base + 1
            } else {
                base
            }
        };

        // Create n_particles neutrons from this product
        for _ in 0..n_particles {
            // Sample outgoing energy and scattering cosine from product distribution
            let (mut e_out, mut mu) = neutron_product.sample(e_in, rng);

            // If scattering is in center-of-mass frame, convert to LAB frame
            // OpenMC: physics.cpp inelastic_scatter() lines 1131-1141
            if reaction.scatter_in_cm {
                let e_cm = e_out;
                let a = awr;

                // Determine outgoing energy in lab frame
                // E = E_cm + (E_in + 2*mu*(A+1)*sqrt(E_in*E_cm)) / ((A+1)^2)
                e_out = e_cm + (e_in + 2.0 * mu * (a + 1.0) * (e_in * e_cm).sqrt())
                        / ((a + 1.0) * (a + 1.0));

                // Determine outgoing angle in lab frame
                // mu = mu * sqrt(E_cm/E) + 1/(A+1) * sqrt(E_in/E)
                mu = mu * (e_cm / e_out).sqrt() + 1.0 / (a + 1.0) * (e_in / e_out).sqrt();

                // Clamp mu to [-1, 1] due to floating point roundoff
                mu = mu.max(-1.0).min(1.0);
            }

            // Create new neutron particle
            let mut new_particle = particle.clone();
            new_particle.energy = e_out;
            crate::inelastic::rotate_direction(&mut new_particle.direction, mu, rng);

            outgoing_neutrons.push(new_particle);
        }
    }

    outgoing_neutrons
}
