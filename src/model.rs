use crate::particle::Particle;
use crate::surface::BoundaryType;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use crate::physics::elastic_scatter;
use crate::data::ATOMIC_WEIGHT_RATIO;
use crate::tally::{Tally, create_tallies_from_specs};
use rayon::prelude::*;
use std::sync::atomic::Ordering;
impl Model {
    pub fn run(&mut self) {
        // println!("Starting particle transport simulation...");

        // Ensure all nuclear data is loaded before transport
        for cell in &mut self.geometry.cells {
            if let Some(material_arc) = &cell.material {
                // Setup phase: clone the inner Material, mutate, and re-wrap in Arc
                let mut material = (**material_arc).clone();
                let _ = material.ensure_nuclides_loaded();
                material.calculate_macroscopic_xs(&vec![1], true);
                cell.material = Some(Arc::new(material));
            }
        }

        // Initialize batch data for all user tallies (updates them in place via RefCell)
        for user_tally_arc in &self.tallies {
            user_tally_arc.initialize_batches(self.settings.batches as usize);
        }

        // Work with the tallies directly during transport
        let tallies = &self.tallies;

        // Get base seed for reproducibility
        let base_seed = self.settings.get_seed();

        for batch in 0..self.settings.batches {
            // println!("Batch {}", batch + 1);

            // Parallelize particle transport within each batch
            // Each particle gets a unique, deterministic seed based on batch and particle index
            (0..self.settings.particles).into_par_iter().for_each(|particle_idx| {
                // Create a unique seed for this particle: base_seed + batch * particles + particle_idx
                let particle_seed = base_seed
                    .wrapping_add((batch as u64).wrapping_mul(self.settings.particles as u64))
                    .wrapping_add(particle_idx as u64);
                let mut rng = StdRng::seed_from_u64(particle_seed);

                // Sample a particle from the source
                let mut particle = self.settings.source.sample(&mut rng);
                particle.alive = true;

                // Transport loop
                while particle.alive {
                    let cell_opt = self.geometry.find_cell((
                        particle.position[0],
                        particle.position[1],
                        particle.position[2],
                    ));
                    if cell_opt.is_none() {
                        panic!("Particle location not found with in any cells at x={}, y={}, z={} - geometry definition error", 
                               particle.position[0], particle.position[1], particle.position[2]);
                    }
                    let cell = cell_opt.unwrap();

                    // Check if cell has material or is void
                    let dist_collision = match &cell.material {
                        Some(material_arc) => {
                            let material = material_arc.as_ref();
                            material
                                .sample_distance_to_collision(particle.energy, &mut rng)
                                .unwrap_or(f64::INFINITY)
                        }
                        None => {
                            // Void cell - no collisions possible
                            // println!("Particle streaming through void cell {:?}", cell.cell_id);
                            f64::INFINITY
                        }
                    };

                    // println!("Sampled distance to collision: {}", dist_collision);

                    // Find closest surface and distance (computed together to avoid recalculation)
                    if let Some((surface_arc, dist_surface)) =
                        cell.closest_surface(particle.position, particle.direction)
                    {
                        // println!("Distance to surface: {}", dist_surface);
                        
                        if dist_surface < dist_collision {
                            // Check boundary type before moving
                            if surface_arc.boundary_type == BoundaryType::Vacuum {
                                // println!(
                                //     "Particle will leak from geometry - killing at {:?}",
                                //     particle.position
                                // );
                                particle.alive = false;
                            } else {
                                // Move slightly past surface to avoid geometric ambiguity
                                const SURFACE_TOLERANCE: f64 = 1e-8;
                                particle.move_by(dist_surface + SURFACE_TOLERANCE);
                                // println!(
                                //     "Particle crossed surface to {:?} (moved {} + tolerance)",
                                //     particle.position, dist_surface
                                // );
                            }
                        } else {
                            // Move to collision point
                            particle.move_by(dist_collision);

                            // We know we have material here because void cells always hit surfaces first
                            let material = cell.material.as_ref().unwrap().as_ref();
                            // Store material_id for filter checking
                            let material_id = material.material_id;
                            // Sample nuclide and reaction
                            let nuclide_name =
                                material.sample_interacting_nuclide(particle.energy, &mut rng);
                            if let Some(nuclide) = material.nuclide_data.get(&nuclide_name) {
                                let reaction = nuclide.sample_reaction(
                                    particle.energy,
                                    &material.temperature,
                                    &mut rng,
                                );
                                if let Some(reaction) = reaction {
                                // println!("Particle collided in cell {:?} at {:?} with nuclide {} via MT {}", cell.cell_id, particle.position, nuclide_name, reaction.mt_number);
                                
                                let mut reaction = reaction;
                                match reaction.mt_number {
                                    2 => {
                                        // Elastic scattering
                                        let awr = *ATOMIC_WEIGHT_RATIO
                                            .get(nuclide_name.as_str())
                                            .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));
                                        elastic_scatter(&mut particle, awr, &mut rng);  // updates particles direction and energy
                                        // println!("Particle elastically scattered at {:?}", particle.position);
                                        // Continue transport (particle.alive remains true)
                                    }
                                    18 => {
                                        // Fission
                                        // println!("Particle caused fission at {:?}", particle.position);
                                        // println!("particle killed as code not ready fission reaction");
                                        // TODO: Sample number of fission neutrons and add to particle bank
                                        particle.alive = false; // Kill original particle for now
                                    }
                                    101 => {
                                        // Absorption (capture)
                                        // Sample which sub-reaction (MT) actually occurred for absorption
                                        let constituent_reaction = nuclide.sample_absorption_constituent(
                                            particle.energy,
                                            &material.temperature,
                                            &mut rng,
                                        );
                                        reaction = constituent_reaction;
                                        particle.alive = false;
                                    }
                                    4 => {
                                        // Inelastic scattering - sample specific constituent reaction
                                        let constituent_reaction = nuclide.sample_inelastic_constituent(
                                            particle.energy,
                                            &material.temperature,
                                            &mut rng,
                                        );
                                        // Overwrite reaction with constituent for tally scoring
                                        reaction = constituent_reaction;
                                        // println!("Particle underwent inelastic scattering at {:?} via MT {} (constituent of MT 4)", 
                                        //         particle.position, constituent_reaction.mt_number);
                                    }
                                    _ => {
                                        // Unknown reaction type - should never happen
                                        panic!("Unknown reaction MT={} at {:?} - sample_reaction returned unexpected MT number", reaction.mt_number, particle.position);
                                    }
                                }
                                
                                // Score user tallies for this reaction (after physics is processed)
                                for tally in tallies.iter() {
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
                        // No surface found - serious geometry error
                        panic!("No surface found for particle at x={}, y={}, z={} with direction [{}, {}, {}] in cell {:?} - geometry definition error",
                            particle.position[0], particle.position[1], particle.position[2],
                            particle.direction[0], particle.direction[1], particle.direction[2],
                            cell.cell_id);
                    }
                }
            }); // End of parallel for_each

            // Update statistics for all user tallies (batch_data already populated via score_event)
            for tally in tallies.iter() {
                tally.update_statistics(self.settings.particles as u32);
            }
        }

        println!("Simulation complete.");
        // Tallies updated in place - no return needed
    }
}
use crate::geometry::Geometry;
// use crate::materials::Materials;
use crate::settings::Settings;
use crate::source::IndependentSource;
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
pub struct Model {
    pub geometry: Geometry,
    pub settings: Settings,
    pub tallies: Vec<Arc<crate::tally::Tally>>,  // Arc allows sharing tallies for parallel updates
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
        let material_arc = Arc::new(material.clone());
        let cell = Cell::new(Some(1), region, Some("sphere_cell".to_string()), Some(material_arc.clone()));
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
        absorption_tally.score = 101; // MT 101 = absorption
        absorption_tally.name = Some("Absorption Tally".to_string());
        absorption_tally.units = "events".to_string();

        let absorption_tally_arc = Arc::new(absorption_tally);

        let mut model = Model {
            geometry,
            settings,
            tallies: vec![Arc::clone(&absorption_tally_arc)]
        };
        assert_eq!(model.settings.particles, 100);
        assert_eq!(model.settings.source.energy, 1e6);
        // Check geometry and material
        assert_eq!(model.geometry.cells.len(), 1);
        let cell_material = &model.geometry.cells[0].material;
        assert!(cell_material.is_some());
        assert!(cell_material
            .as_ref()
            .unwrap()
            .nuclides
            .contains_key("Li6"));
        // Run the model and ensure it executes without panicking
        model.run();

        // Verify the absorption tally was updated in place
        assert_eq!(absorption_tally_arc.name, Some("Absorption Tally".to_string()));
        assert_eq!(absorption_tally_arc.units, "events");
        assert_eq!(absorption_tally_arc.n_batches.load(Ordering::Relaxed), 10);

        println!("Test tally results:");
        println!("Absorption Tally: {}", absorption_tally_arc);

        println!("✓ Absorption tally verified successfully - tally updated in place!");
    }
}
