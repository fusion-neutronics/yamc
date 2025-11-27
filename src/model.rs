use crate::data::ATOMIC_WEIGHT_RATIO;
use crate::geometry::Geometry;
use crate::inelastic::inelastic_scatter;
use crate::particle::Particle;
use crate::physics::elastic_scatter;
use crate::reaction::Reaction;
use crate::settings::Settings;
use crate::source::IndependentSource;
use crate::surface::BoundaryType;
use crate::tally::{create_tallies_from_specs, Tally};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rand_pcg::Pcg64;
use rayon::prelude::*;
use std::collections::VecDeque;
use std::sync::{Arc, atomic::Ordering};

#[derive(Debug, Clone)]
pub struct Model {
    pub geometry: Geometry,
    pub settings: Settings,
    pub tallies: Vec<Arc<crate::tally::Tally>>, // Arc allows sharing tallies for parallel updates
}

impl Model {
    /// Create a new Model with particle banking system
    pub fn new(geometry: Geometry, settings: Settings, tallies: Vec<Arc<crate::tally::Tally>>) -> Self {
        Model {
            geometry,
            settings,
            tallies,
        }
    }

    pub fn run(&mut self) {
        // Ensure all nuclear data is loaded before transport
        for cell in &mut self.geometry.cells {
            if let Some(material_arc) = &cell.material {
                let mut material = (**material_arc).clone();
                let _ = material.ensure_nuclides_loaded();
                material.calculate_macroscopic_xs(&vec![1], true);
                cell.material = Some(Arc::new(material));
            }
        }

        // Ensure all tallies are initialized with the correct number of batches
        for tally_arc in &self.tallies {
            tally_arc.initialize_batches_shared(self.settings.batches as usize);
        }

        let tallies = &self.tallies;
        let base_seed = self.settings.get_seed();

        for batch in 0..self.settings.batches {
            println!("Starting batch {}/{}", batch + 1, self.settings.batches);
            // Parallel execution per batch with particle banking for multi-neutron reactions
            (0..self.settings.particles).into_par_iter().for_each(|particle_idx| {
                // OpenMC-style per-particle seeding with stride
                // Each particle gets a unique, reproducible seed
                const PARTICLE_STRIDE: u64 = 152917; // Same stride as OpenMC for compatibility
                let particle_seed = base_seed
                    .wrapping_add((batch as u64).wrapping_mul(self.settings.particles as u64).wrapping_mul(PARTICLE_STRIDE))
                    .wrapping_add((particle_idx as u64).wrapping_mul(PARTICLE_STRIDE));
                
                // Use PCG64 (fast, high-quality RNG similar to OpenMC's PCG)
                let mut rng = Pcg64::seed_from_u64(particle_seed);

                // Local particle queue for this thread (for secondary particles)
                let mut particle_queue: VecDeque<Particle> = VecDeque::new();
                
                // Start with the source particle
                let mut particle = self.settings.source.sample(&mut rng);
                particle.alive = true;
                particle_queue.push_back(particle);
                
                // Process all particles in the queue (primary + secondaries)
                let mut particles_processed = 0;
                const MAX_PARTICLES_PER_HISTORY: usize = 1000; // Safety limit
                
                while let Some(mut particle) = particle_queue.pop_front() {
                    particles_processed += 1;
                    if particles_processed > MAX_PARTICLES_PER_HISTORY {
                        break; // Prevent infinite particle multiplication
                    }

                while particle.alive {
                    // Use cached cell index if available, otherwise search for it
                    let cell_index = match particle.current_cell_index {
                        Some(idx) => idx,
                        None => {
                            // Need to find the cell (after initialization or surface crossing)
                            match self.geometry.find_cell_index((
                                particle.position[0],
                                particle.position[1],
                                particle.position[2],
                            )) {
                                Some(idx) => {
                                    particle.current_cell_index = Some(idx);
                                    idx
                                }
                                None => {
                                    panic!(
                                        "Particle location not found within any cells at x={}, y={}, z={} - geometry definition error",
                                        particle.position[0], particle.position[1], particle.position[2]
                                    );
                                }
                            }
                        }
                    };
                    let cell = &self.geometry.cells[cell_index];

                    let dist_collision = match &cell.material {
                        Some(material_arc) => {
                            let material = material_arc.as_ref();
                            material
                                .sample_distance_to_collision(particle.energy, &mut rng)
                                .unwrap_or(f64::INFINITY)
                        }
                        None => {
                            f64::INFINITY
                        }
                    };

                    if let Some((surface_arc, dist_surface)) =
                        cell.closest_surface(particle.position, particle.direction)
                    {
                        if dist_surface < dist_collision {
                            if surface_arc.boundary_type == BoundaryType::Vacuum {
                                particle.alive = false;
                            } else {
                                const SURFACE_TOLERANCE: f64 = 1e-8;
                                particle.move_by(dist_surface + SURFACE_TOLERANCE);
                                // Invalidate cached cell after surface crossing
                                particle.current_cell_index = None;
                            }
                        } else {
                            particle.move_by(dist_collision);
                            let material = cell.material.as_ref().unwrap().as_ref();
                            let material_id = material.material_id;
                            let nuclide_name =
                                material.sample_interacting_nuclide(particle.energy, &mut rng);
                            if let Some(nuclide) = material.nuclide_data.get(&nuclide_name) {
                                let reaction = nuclide.sample_reaction_no_autoload(
                                    particle.energy,
                                    &material.temperature,
                                    &mut rng,
                                );
                                if let Some(reaction) = reaction {
                                    let mut reaction = reaction;
                                    match reaction.mt_number {
                                        1001 => {
                                            // Scattering - sample which constituent MT actually occurred
                                            let constituent_reaction = nuclide.sample_scattering_constituent(
                                                particle.energy,
                                                &material.temperature,
                                                &mut rng,
                                            );

                                            // Handle the constituent reaction
                                            match constituent_reaction.mt_number {
                                                2 => {
                                                    // Elastic scatter - no products to handle
                                                    let awr = *ATOMIC_WEIGHT_RATIO
                                                        .get(nuclide_name.as_str())
                                                        .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));
                                                    elastic_scatter(&mut particle, awr, &mut rng);
                                                }
                                                50..=91 => {
                                                    // Inelastic scatter (MT 50-91) - always produces exactly 1 neutron
                                                    let outgoing_particles = crate::inelastic::analytical_inelastic_scatter(
                                                        &particle, &constituent_reaction, &nuclide_name, &mut rng
                                                    );

                                                    // MT 50-91 should always return exactly 1 particle
                                                    assert_eq!(
                                                        outgoing_particles.len(),
                                                        1,
                                                        "MT {} should produce exactly 1 neutron, got {}",
                                                        constituent_reaction.mt_number,
                                                        outgoing_particles.len()
                                                    );

                                                    let outgoing = &outgoing_particles[0];
                                                    particle.energy = outgoing.energy;
                                                    particle.direction = outgoing.direction;
                                                    particle.position = outgoing.position;
                                                }
                                                _ => {
                                                    // Other scattering reactions - handle products, multiplicity, banking
                                                    // Use scatter module for general scattering (n,2n), (n,3n), (n,n'alpha), etc.
                                                    let outgoing_particles = crate::scatter::scatter(
                                                        &particle, &constituent_reaction, &nuclide_name, &mut rng
                                                    );

                                                    // All scattering MTs should produce at least 1 neutron
                                                    assert!(
                                                        !outgoing_particles.is_empty(),
                                                        "Scattering MT {} produced no neutrons - this should not happen",
                                                        constituent_reaction.mt_number
                                                    );

                                                    if outgoing_particles.len() == 1 {
                                                        // Single neutron - update current particle
                                                        let outgoing = &outgoing_particles[0];
                                                        particle.energy = outgoing.energy;
                                                        particle.direction = outgoing.direction;
                                                        particle.position = outgoing.position;
                                                    } else {
                                                        // Multi-neutron reaction - bank secondary particles
                                                        for (i, outgoing_particle) in outgoing_particles.into_iter().enumerate() {
                                                            if i == 0 {
                                                                particle.energy = outgoing_particle.energy;
                                                                particle.direction = outgoing_particle.direction;
                                                                particle.position = outgoing_particle.position;
                                                            } else {
                                                                particle_queue.push_back(outgoing_particle);
                                                            }
                                                        }
                                                    }
                                                }
                                            }

                                            // Use constituent reaction for tally scoring
                                            reaction = constituent_reaction;
                                        }
                                        18 => {
                                            particle.alive = false;
                                        }
                                        101 => {
                                            let constituent_reaction = nuclide.sample_absorption_constituent(
                                                particle.energy,
                                                &material.temperature,
                                                &mut rng,
                                            );
                                            reaction = constituent_reaction;
                                            particle.alive = false;
                                        }
                                        _ => {
                                            // For unknown reactions, just absorb the particle
                                            particle.alive = false;
                                        }
                                    }                                    for tally in tallies.iter() {
                                        tally.score_event(reaction.mt_number, cell, material_id, batch as usize);
                                    }
                                } else {
                                    panic!(
                                        "No valid reaction found for nuclide {} at energy {}",
                                        nuclide_name, particle.energy
                                    );
                                }
                            } else {
                                panic!("Nuclide {} not found in material data", nuclide_name);
                            }
                        }
                    } else {
                        panic!("No surface found for particle at x={}, y={}, z={} with direction [{}, {}, {}] in cell {:?} - geometry definition error",
                            particle.position[0], particle.position[1], particle.position[2],
                            particle.direction[0], particle.direction[1], particle.direction[2],
                            cell.cell_id);
                    }
                } // end of particle transport while loop
                } // end of particle queue processing while loop
            }); // end of parallel for_each

            for tally in tallies.iter() {
                tally.update_statistics(self.settings.particles as u32);
            }
        }

        println!("Simulation complete.");
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;
    use crate::geometry::Geometry;
    use crate::material::Material;
    use crate::region::{HalfspaceType, Region};
    use crate::settings::Settings;
    use crate::source::IndependentSource;
    use crate::surface::{BoundaryType, Surface, SurfaceKind};
    // No duplicate import

    #[test]
    fn test_model_construction() {
        // Sphere surface with vacuum boundary
        let sphere = Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 2.0,
            },
            boundary_type: BoundaryType::Vacuum,
        };
        let region = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere)));
        // Material with Li6 nuclide
        let mut material = Material::new();
        material.set_density("g/cc", 1.0).unwrap(); // Add density to fix test
        material.add_nuclide("Li6", 1.0).unwrap();
        let mut nuclide_json_map = std::collections::HashMap::new();
        nuclide_json_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
        material.read_nuclides_from_json(&nuclide_json_map).unwrap();
        // Prepare material cross sections before creating the Arc
        material.calculate_macroscopic_xs(&vec![1], true);
        let material_arc = Arc::new(material.clone());
        let cell = Cell::new(
            Some(1),
            region,
            Some("sphere_cell".to_string()),
            Some(material_arc.clone()),
        );
        let geometry = Geometry { cells: vec![cell] };
        let source = IndependentSource {
            space: [0.0, 0.0, 0.0],
            angle: crate::stats::AngularDistribution::new_monodirectional(0.0, 0.0, 1.0),
            energy: 1e6,
        };
        let settings = Settings {
            particles: 100,
            batches: 10,
            source: source.clone(),
            seed: Some(1),
        };

        // Create a tally for absorption reactions with a custom name
        let mut absorption_tally = crate::tally::Tally::new();
        absorption_tally.scores = vec![101]; // MT 101 = absorption
        absorption_tally.name = Some("Absorption Tally".to_string());
        absorption_tally.units = "events".to_string();
        absorption_tally.initialize_batches(settings.batches as usize);
        let absorption_tally_arc = Arc::new(absorption_tally);

        let mut model = Model::new(
            geometry,
            settings,
            vec![Arc::clone(&absorption_tally_arc)],
        );
        assert_eq!(model.settings.particles, 100);
        assert_eq!(model.settings.source.energy, 1e6);
        // Check geometry and material
        assert_eq!(model.geometry.cells.len(), 1);
        let cell_material = &model.geometry.cells[0].material;
        assert!(cell_material.is_some());
        assert!(cell_material.as_ref().unwrap().nuclides.contains_key("Li6"));
        // Run the model and ensure it executes without panicking
        model.run();

        // Verify the absorption tally was updated in place
        assert_eq!(
            absorption_tally_arc.name,
            Some("Absorption Tally".to_string())
        );
        assert_eq!(absorption_tally_arc.units, "events");
        assert_eq!(absorption_tally_arc.n_batches.load(Ordering::Relaxed), 10);

        println!("Test tally results:");
        println!("Absorption Tally: {}", absorption_tally_arc);

        println!("✓ Absorption tally verified successfully - tally updated in place!");
    }
}
