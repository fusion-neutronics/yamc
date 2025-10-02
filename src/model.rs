use crate::particle::Particle;
use crate::surface::BoundaryType;
use rand::Rng;
use crate::physics::elastic_scatter;
use crate::data::ATOMIC_WEIGHT_RATIO;
impl Model {
    pub fn run(&self) {
        println!("Starting particle transport simulation...");

        // Ensure all nuclear data is loaded before transport
        for cell in &self.geometry.cells {
            if let Some(material_arc) = &cell.material {
                let mut material = material_arc.lock().unwrap();
                let _ = material.ensure_nuclides_loaded();
                material.calculate_macroscopic_xs(&vec![1], true);
            }
        }

        let mut rng = rand::thread_rng();
        for batch in 0..self.settings.batches {
            println!("Batch {}", batch + 1);
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
                        println!("Particle leaked from geometry at {:?}", particle.position);
                        particle.alive = false;
                        break;
                    }
                    let cell = cell_opt.unwrap();

                    // Get material for this cell
                    let material = match &cell.material {
                        Some(mat_arc_mutex) => {
                            let mat = mat_arc_mutex.lock().unwrap();
                            mat
                        }
                        None => {
                            println!("No material in cell {}", cell.cell_id);
                            particle.alive = false;
                            break;
                        }
                    };

                    // Sample distance to collision
                    let dist_collision = material
                        .sample_distance_to_collision(particle.energy, &mut rng)
                        .unwrap_or(f64::INFINITY);
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
                            .unwrap_or(f64::INFINITY);
                        if dist_surface < dist_collision {
                            // Move to surface
                            for i in 0..3 {
                                particle.position[i] += particle.direction[i] * dist_surface;
                            }
                            // Check boundary type
                            if surface_arc.boundary_type == BoundaryType::Vacuum {
                                println!(
                                    "Particle leaked from geometry at {:?}",
                                    particle.position
                                );
                                particle.alive = false;
                            } else {
                                // For now, just kill the particle at any boundary
                                println!(
                                    "Particle hit non-vacuum boundary at {:?}",
                                    particle.position
                                );
                                particle.alive = false;
                            }
                        } else {
                            // Move to collision point
                            for i in 0..3 {
                                particle.position[i] += particle.direction[i] * dist_collision;
                            }
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
                                    println!("Particle collided in cell {} at {:?} with nuclide {} via MT {}", cell.cell_id, particle.position, nuclide_name, reaction.mt_number);
                                    
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
                                            // TODO: Sample number of fission neutrons and add to particle bank
                                            particle.alive = false; // Kill original particle for now
                                        }
                                        101 => {
                                            // Absorption (capture)
                                            println!("Particle absorbed at {:?} (MT=101 absorption)", particle.position);
                                            particle.alive = false;
                                        }
                                        3 => {
                                            // Nonelastic scattering (inelastic + other)
                                            println!("Particle underwent nonelastic scattering at {:?}", particle.position);
                                            // TODO: Sample outgoing energy and angle from nuclear data
                                            particle.alive = false; // Kill particle for now until inelastic implemented
                                        }
                                        _ => {
                                            // Unknown reaction type - should never happen
                                            panic!("Unknown reaction MT={} at {:?} - sample_reaction returned unexpected MT number", reaction.mt_number, particle.position);
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
                            particle.alive = false; // End after one collision for now
                        }
                    } else {
                        // No surface found, particle leaks
                        println!("Particle leaked from geometry at {:?}", particle.position);
                        particle.alive = false;
                    }
                }
            }
        }
        println!("Simulation complete.");
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
            surface_id: 1,
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
        let cell = Cell {
            cell_id: 1,
            name: Some("sphere_cell".to_string()),
            region,
            material: Some(material_arc.clone()),
        };
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
        let model = Model { geometry, settings };
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
        model.run();
    }
}
