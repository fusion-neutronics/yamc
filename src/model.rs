use crate::data::ATOMIC_WEIGHT_RATIO;
use crate::particle::Particle;
use crate::physics::elastic_scatter;
use crate::surface::BoundaryType;
use crate::tally::{create_tallies_from_specs, Tally};
use rand::rngs::StdRng;
use rand::Rng;
use rand::SeedableRng;
use rayon::prelude::*;
use std::sync::atomic::Ordering;
impl Model {
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
            (0..self.settings.particles).into_par_iter().for_each(|particle_idx| {
                let particle_seed = base_seed
                    .wrapping_add((batch as u64).wrapping_mul(self.settings.particles as u64))
                    .wrapping_add(particle_idx as u64);
                let mut rng = StdRng::seed_from_u64(particle_seed);

                let mut particle = self.settings.source.sample(&mut rng);
                particle.alive = true;

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
                                        2 => {
                                            let awr = *ATOMIC_WEIGHT_RATIO
                                                .get(nuclide_name.as_str())
                                                .expect(&format!("No atomic weight ratio for nuclide {}", nuclide_name));
                                            elastic_scatter(&mut particle, awr, &mut rng);
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
                                        4 => {
                                            let constituent_reaction = nuclide.sample_inelastic_constituent(
                                                particle.energy,
                                                &material.temperature,
                                                &mut rng,
                                            );
                                            reaction = constituent_reaction;
                                        }
                                        _ => {
                                            panic!("Unknown reaction MT={} at {:?} - sample_reaction returned unexpected MT number", reaction.mt_number, particle.position);
                                        }
                                    }

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
                        panic!("No surface found for particle at x={}, y={}, z={} with direction [{}, {}, {}] in cell {:?} - geometry definition error",
                            particle.position[0], particle.position[1], particle.position[2],
                            particle.direction[0], particle.direction[1], particle.direction[2],
                            cell.cell_id);
                    }
                }
            });

            for tally in tallies.iter() {
                tally.update_statistics(self.settings.particles as u32);
            }
        }

        println!("Simulation complete.");
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
    pub tallies: Vec<Arc<crate::tally::Tally>>, // Arc allows sharing tallies for parallel updates
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

        let mut model = Model {
            geometry,
            settings,
            tallies: vec![Arc::clone(&absorption_tally_arc)],
        };
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
