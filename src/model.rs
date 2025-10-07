use crate::particle::Particle;
use crate::surface::BoundaryType;
use rand::Rng;
use crate::physics::elastic_scatter;
use crate::data::ATOMIC_WEIGHT_RATIO;
use crate::tally::{Tally, create_tallies_from_specs};
impl Model {
    /// Create a new Model with pre-computed optimization flags
    pub fn new(
        geometry: Geometry,
        settings: Settings,
        tallies: Vec<crate::tally::Tally>,
    ) -> Self {
        // Check if any tallies require detailed absorption sampling
        let needs_detailed_absorption = tallies.iter().any(|t| {
            // MT 101 is absorption, but we need detailed sampling for its constituents
            let absorption_constituents = [
                102, 103, 104, 105, 106, 107, 108, 109, 111, 112, 113, 114, 115, 116, 117,
                155, 182, 191, 192, 193, 197,
            ];
            absorption_constituents.contains(&t.score)
        });

        // Check if any tallies require detailed nonelastic sampling
        let needs_detailed_nonelastic = tallies.iter().any(|t| {
            // MT 3 is nonelastic, but we need detailed sampling for its constituents
            let nonelastic_constituents = [
                4, 5, 11, 16, 17, 22, 23, 24, 25, 27, 28, 29, 30, 32, 33, 34, 35,
                36, 37, 41, 42, 44, 45, 152, 153, 154, 156, 157, 158, 159, 160,
                161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172,
                173, 174, 175, 176, 177, 178, 179, 180, 181, 183, 184, 185,
                186, 187, 188, 189, 190, 194, 195, 196, 198, 199, 200,
            ];
            nonelastic_constituents.contains(&t.score)
        });

        Model {
            geometry,
            settings,
            tallies,
            needs_detailed_absorption,
            needs_detailed_nonelastic,
        }
    }

    pub fn run(&self) -> Vec<Tally> {
        println!("Starting particle transport simulation...");

        // Ensure all nuclear data is loaded before transport
        for cell in &self.geometry.cells {
            if let Some(material_arc) = &cell.material {
                let mut material = material_arc.lock().unwrap();
                let _ = material.ensure_nuclides_loaded();
                material.calculate_macroscopic_xs(&vec![1], true);
            }
        }

        // Initialize tallies from user specifications  
        let mut tallies = create_tallies_from_specs(&self.tallies);

        let mut rng = rand::thread_rng();
        for batch in 0..self.settings.batches {
            println!("Batch {}", batch + 1);
            
            // Initialize leakage counter for this batch
            let mut batch_leakage = 0u32;
            
            // Initialize user tally counters for this batch (excluding leakage)
            let mut user_batch_counts: Vec<u32> = vec![0; tallies.len() - 1];
            
            for _ in 0..self.settings.particles {
                // Sample a particle from the source via settings
                let mut particle = self.settings.source.sample();
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
                        Some(mat_arc_mutex) => {
                            let material = mat_arc_mutex.lock().unwrap();
                            material
                                .sample_distance_to_collision(particle.energy, &mut rng)
                                .unwrap_or(f64::INFINITY)
                        }
                        None => {
                            // Void cell - no collisions possible
                            println!("Particle streaming through void cell {:?}", cell.cell_id);
                            f64::INFINITY
                        }
                    };

                    println!("Sampled distance to collision: {}", dist_collision);

                    // Find closest surface and distance
                    if let Some(surface_arc) =
                        cell.closest_surface(particle.position, particle.direction)
                    {
                        let dist_surface = surface_arc
                            .distance_to_surface(
                                [
                                    particle.position[0],
                                    particle.position[1],
                                    particle.position[2],
                                ],
                                particle.direction,
                            )
                            .unwrap_or_else(|| {
                                panic!("Failed to calculate distance to surface for particle at x={}, y={}, z={} with direction [{}, {}, {}] in cell {:?}",
                                    particle.position[0], particle.position[1], particle.position[2],
                                    particle.direction[0], particle.direction[1], particle.direction[2],
                                    cell.cell_id);
                            });
                        
                        println!("Distance to surface: {}", dist_surface);
                        
                        if dist_surface < dist_collision {
                            // Check boundary type before moving
                            if surface_arc.boundary_type == BoundaryType::Vacuum {
                                println!(
                                    "Particle will leak from geometry - killing at {:?}",
                                    particle.position
                                );
                                particle.alive = false;
                                batch_leakage += 1; // Count this leakage
                            } else {
                                // Move slightly past surface to avoid geometric ambiguity
                                const SURFACE_TOLERANCE: f64 = 1e-8;
                                particle.move_by(dist_surface + SURFACE_TOLERANCE);
                                println!(
                                    "Particle crossed surface to {:?} (moved {} + tolerance)",
                                    particle.position, dist_surface
                                );
                            }
                        } else {
                            // Move to collision point
                            particle.move_by(dist_collision);

                            // We know we have material here because void cells always hit surfaces first
                            let material = cell.material.as_ref().unwrap().lock().unwrap();
                            // Store material_id for filter checking (to avoid double-locking later)
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
                                println!("Particle collided in cell {:?} at {:?} with nuclide {} via MT {}", cell.cell_id, particle.position, nuclide_name, reaction.mt_number);
                                
                                // Determine the specific MT for tallying (use detailed sampling if needed)
                                let tally_mt = match reaction.mt_number {
                                    101 if self.needs_detailed_absorption => {
                                        // Sample specific absorption subreaction for tallying
                                        material.sample_absorption_subreaction(&nuclide_name, particle.energy, &mut rng)
                                            .unwrap_or(101) // Fall back to MT 101 if no specific subreaction
                                    }
                                    3 if self.needs_detailed_nonelastic => {
                                        // Sample specific nonelastic subreaction for tallying
                                        material.sample_nonelastic_subreaction(&nuclide_name, particle.energy, &mut rng)
                                            .unwrap_or(3) // Fall back to MT 3 if no specific subreaction
                                    }
                                    _ => reaction.mt_number, // Use the sampled MT as-is
                                };

                                match reaction.mt_number {
                                    2 => {
                                        // Elastic scattering
                                        let awr = *ATOMIC_WEIGHT_RATIO
                                            .get(nuclide_name.as_str())
                                            .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));
                                        elastic_scatter(&mut particle, awr, &mut rng);  // updates particles direction and energy
                                        println!("Particle elastically scattered at {:?}", particle.position);
                                        // Continue transport (particle.alive remains true)
                                    }
                                    18 => {
                                        // Fission
                                        println!("Particle caused fission at {:?}", particle.position);
                                        println!("particle killed as code not ready fission reaction");
                                        // TODO: Sample number of fission neutrons and add to particle bank
                                        particle.alive = false; // Kill original particle for now
                                    }
                                    101 => {
                                        // Absorption (capture)
                                        println!("Particle absorbed at {:?} (MT={} specific subreaction)", particle.position, tally_mt);
                                        particle.alive = false;
                                    }
                                    3 => {
                                        // Nonelastic scattering (inelastic + other)
                                        println!("Particle underwent nonelastic scattering at {:?} (MT={} specific subreaction)", particle.position, tally_mt);
                                        println!("particle killed as code not ready nonelastic reaction");
                                        // TODO: Sample outgoing energy and angle from nuclear data
                                        particle.alive = false; // Kill particle for now until inelastic implemented
                                    }
                                    _ => {
                                        // Unknown reaction type - should never happen
                                        panic!("Unknown reaction MT={} at {:?} - sample_reaction returned unexpected MT number", reaction.mt_number, particle.position);
                                    }
                                }
                                
                                // Score user tallies for this reaction (after physics is processed)
                                for (i, tally_spec) in self.tallies.iter().enumerate() {
                                    if tally_spec.score == tally_mt {
                                        // Check if this event passes all filters for this tally
                                        let passes_filters = if tally_spec.filters.is_empty() {
                                            // No filters means score all events
                                            true
                                        } else {
                                            // Check if ALL filters match this event (intersection of filters)
                                            tally_spec.filters.iter().all(|filter| match filter {
                                                crate::tally::Filter::Cell(cell_filter) => {
                                                    cell.cell_id.map_or(false, |id| cell_filter.matches(id))
                                                },
                                                crate::tally::Filter::Material(material_filter) => {
                                                    // Check if the cell's material matches this material filter
                                                    // Use the stored material_id to avoid double-locking the mutex
                                                    material_filter.matches(material_id)
                                                }
                                            })
                                        };
                                        
                                        if passes_filters {
                                            user_batch_counts[i] += 1;
                                        }
                                    }
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
            }
            
            // Store batch leakage at end of each batch (first tally is always leakage)
            tallies[0].add_batch(batch_leakage, self.settings.particles as u32);
            
            // Store user tally results for this batch (starting from index 1)
            for (i, count) in user_batch_counts.iter().enumerate() {
                tallies[i + 1].add_batch(*count, self.settings.particles as u32);
            }
        }
        
        println!("Simulation complete.");
        tallies
    }
}
use crate::geometry::Geometry;
// use crate::materials::Materials;
use crate::settings::Settings;
use crate::source::Source;
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
pub struct Model {
    pub geometry: Geometry,
    pub settings: Settings,
    pub tallies: Vec<crate::tally::Tally>,
    // Pre-computed flags for performance optimization
    needs_detailed_absorption: bool,
    needs_detailed_nonelastic: bool,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::Cell;
    use crate::geometry::Geometry;
    use crate::material::Material;
    use crate::region::{HalfspaceType, Region};
    use crate::settings::Settings;
    use crate::source::Source;
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
        let material_arc = Arc::new(Mutex::new(material));
        let cell = Cell::new(Some(1), region, Some("sphere_cell".to_string()), Some(material_arc.clone()));
        let geometry = Geometry { cells: vec![cell] };
        let source = Source {
            position: [0.0, 0.0, 0.0],
            direction: [0.0, 0.0, 1.0],
            energy: 1e6,
        };
        let settings = Settings {
            particles: 100,
            batches: 10,
            source: source.clone(),
        };
        
        // Create a tally for absorption reactions with a custom name
        let mut absorption_tally = crate::tally::Tally::new();
        absorption_tally.score = 101; // MT 101 = absorption
        absorption_tally.name = Some("Absorption Tally".to_string());
        
        let model = Model::new(geometry, settings, vec![absorption_tally]);
        assert_eq!(model.settings.particles, 100);
        assert_eq!(model.settings.source.energy, 1e6);
        // Check geometry and material
        assert_eq!(model.geometry.cells.len(), 1);
        let cell_material = &model.geometry.cells[0].material;
        assert!(cell_material.is_some());
        assert!(cell_material
            .as_ref()
            .unwrap()
            .lock()
            .unwrap()
            .nuclides
            .contains_key("Li6"));
        // Run the model and ensure it executes without panicking
        let tallies = model.run();

        // Verify we got the expected number of tallies (leakage + 1 user tally)
        assert_eq!(tallies.len(), 2, "Expected 2 tallies: leakage + absorption");
        
        // Verify leakage tally (index 0)
        assert_eq!(tallies[0].name, Some("Leakage".to_string()));
        assert_eq!(tallies[0].units, "particles");
        assert_eq!(tallies[0].n_batches, 10);
        
        // Verify absorption tally (index 1) 
        assert_eq!(tallies[1].name, Some("Absorption Tally".to_string()));
        assert_eq!(tallies[1].units, "events");
        assert_eq!(tallies[1].n_batches, 10);
        
        println!("Test tally results:");
        for (i, tally) in tallies.iter().enumerate() {
            println!("Tally {}: {}", i, tally);
        }
        
        println!("✓ Both leakage and absorption tallies verified successfully");
    }
}
