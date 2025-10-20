use crate::particle::Particle;
use crate::surface::BoundaryType;
use rand::Rng;
use crate::physics::elastic_scatter;
use crate::data::ATOMIC_WEIGHT_RATIO;
use crate::tally::{Tally, create_tallies_from_specs};
impl Model {
    pub fn run(&self) -> Vec<Tally> {
        // println!("Starting particle transport simulation...");

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
            // println!("Batch {}", batch + 1);
            
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
                            // println!("Particle streaming through void cell {:?}", cell.cell_id);
                            f64::INFINITY
                        }
                    };

                    // println!("Sampled distance to collision: {}", dist_collision);

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
                        
                        // println!("Distance to surface: {}", dist_surface);
                        
                        if dist_surface < dist_collision {
                            // Check boundary type before moving
                            if surface_arc.boundary_type == BoundaryType::Vacuum {
                                println!("[DEBUG] Particle leaked at position {:?}", particle.position);
                                particle.alive = false;
                                batch_leakage += 1; // Count this leakage
                                break; // Exit transport loop for this particle
                            } else {
                                // Move slightly past surface to avoid geometric ambiguity
                                const SURFACE_TOLERANCE: f64 = 1e-8;
                                particle.move_by(dist_surface + SURFACE_TOLERANCE);
                            }
                                        println!("[DEBUG] Particle absorbed in cell {:?} at position {:?}", cell.cell_id, particle.position);
                                        break; // Exit transport loop for this particle
                                        println!("[DEBUG] Particle absorbed in cell {:?} at position {:?}", cell.cell_id, particle.position);
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
                                        panic!("Fission (MT=18) encountered at {:?} for nuclide {}. Fission is not handled in this model.", particle.position, nuclide_name);
                                    }
                                    101 => {
                                        // Absorption (capture)
                                        // Sample which sub-reaction (MT) actually occurred for absorption
                                        let constituent_reaction = nuclide.sample_absorption_constituent(
                                            particle.energy,
                                            &material.temperature,
                                            &mut rng,
                                        );
                                        // println!(
                                        //     "Particle absorbed at {:?} (MT=101 absorption), sub-reaction MT={}",
                                        //     particle.position,
                                        //     constituent_reaction.mt_number
                                        // );
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
                                // Absorption event logic: MT 101 or any constituent absorption MT
                                let is_absorption_mt = reaction.mt_number == 101
                                    || matches!(reaction.mt_number,
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
                                                    // println!("Particle streaming through void cell {:?}", cell.cell_id);
                                                    f64::INFINITY
                                                }
                                            };

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

                                                if dist_surface < dist_collision {
                                                    if surface_arc.boundary_type == BoundaryType::Vacuum {
                                                        println!("[DEBUG] Particle leaked at position {:?}", particle.position);
                                                        particle.alive = false;
                                                        batch_leakage += 1; // Count this leakage
                                                        break; // Exit transport loop for this particle
                                                    } else {
                                                        // Move slightly past surface to avoid geometric ambiguity
                                                        const SURFACE_TOLERANCE: f64 = 1e-8;
                                                        particle.move_by(dist_surface + SURFACE_TOLERANCE);
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
                                                            let mut reaction = reaction;
                                                            match reaction.mt_number {
                                                                2 => {
                                                                    // Elastic scattering
                                                                    let awr = *ATOMIC_WEIGHT_RATIO
                                                                        .get(nuclide_name.as_str())
                                                                        .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));
                                                                    elastic_scatter(&mut particle, awr, &mut rng);  // updates particles direction and energy
                                                                    // Continue transport (particle.alive remains true)
                                                                }
                                                                18 => {
                                                                    panic!("Fission (MT=18) encountered at {:?} for nuclide {}. Fission is not handled in this model.", particle.position, nuclide_name);
                                                                }
                                                                101 => {
                                                                    // Absorption (capture)
                                                                    let constituent_reaction = nuclide.sample_absorption_constituent(
                                                                        particle.energy,
                                                                        &material.temperature,
                                                                        &mut rng,
                                                                    );
                                                                    reaction = constituent_reaction;
                                                                    particle.alive = false;
                                                                    break; // Exit transport loop for this particle
                                                                }
                                                                4 => {
                                                                    // Inelastic scattering - sample specific constituent reaction
                                                                    let constituent_reaction = nuclide.sample_inelastic_constituent(
                                                                        particle.energy,
                                                                        &material.temperature,
                                                                        &mut rng,
                                                                    );
                                                                    reaction = constituent_reaction;
                                                                    // Continue transport (particle.alive remains true)
                                                                }
                                                                _ => {
                                                                    panic!("Unknown reaction MT={} at {:?} - sample_reaction returned unexpected MT number", reaction.mt_number, particle.position);
                                                                }
                                                            }

                                                            // Score user tallies for this reaction (after physics is processed)
                                                            let is_absorption_mt = reaction.mt_number == 101
                                                                || matches!(reaction.mt_number,
                                                                    102 | 103 | 104 | 105 | 106 | 107 | 108 | 109 | 111 | 112 | 113 | 114 |
                                                                    115 | 116 | 117 | 155 | 182 | 191 | 192 | 193 | 197
                                                                )
                                                                || (600..650).contains(&reaction.mt_number)
                                                                || (650..700).contains(&reaction.mt_number)
                                                                || (700..750).contains(&reaction.mt_number)
                                                                || (750..800).contains(&reaction.mt_number)
                                                                || (800..850).contains(&reaction.mt_number);

                                                            if is_absorption_mt {
                                                                for (j, t) in tallies[1..].iter_mut().enumerate() {
                                                                    if t.score == 101 {
                                                                        if t.score_event(reaction.mt_number, cell, material_id) {
                                                                            user_batch_counts[j] += 1;
                                                                        }
                                                                    }
                                                                }
                                                            } else {
                                                                for (i, tally) in tallies[1..].iter_mut().enumerate() {
                                                                    if tally.score_event(reaction.mt_number, cell, material_id) {
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
                                                panic!("No surface found for particle at x={}, y={}, z={} with direction [{}, {}, {}] in cell {:?} - geometry definition error",
                                                    particle.position[0], particle.position[1], particle.position[2],
                                                    particle.direction[0], particle.direction[1], particle.direction[2],
                                                    cell.cell_id);
                                            }
                                        }
                                        // After transport loop, check that particle is not still alive (i.e., must have been absorbed or leaked)
                                        if particle.alive {
                                            panic!("Particle survived transport loop without being absorbed or leaked. This indicates a logic error. Position: {:?}, Direction: {:?}", particle.position, particle.direction);
                                        }
