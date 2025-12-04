use yaml::*;
use yaml::model::Model;
use yaml::source::IndependentSource;
use yaml::stats::AngularDistribution;
use yaml::settings::Settings;
use yaml::tally::{Score, Tally};
use yaml::filter::Filter;
use yaml::filters::{CellFilter, EnergyFilter};
use std::sync::Arc;
use std::time::Instant;

fn main() {
    // Create a simple spherical geometry with Li6
    let sphere1 = Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 1.0,
        },
        boundary_type: BoundaryType::Transmission,
    };
    
    let sphere2 = Surface {
        surface_id: Some(2),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 200.0,
        },
        boundary_type: BoundaryType::Vacuum,
    };

    let region1 = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere1.clone())));
    
    let sphere_for_region2 = Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 1.0,
        },
        boundary_type: BoundaryType::Transmission,
    };
    
    let region2_inner = Region::new_from_halfspace(HalfspaceType::Above(Arc::new(sphere_for_region2)));
    let region2_outer = Region::new_from_halfspace(HalfspaceType::Below(Arc::new(sphere2)));
    let region2 = region2_inner.intersection(&region2_outer);

    // Create Li6 material
    let mut material = Material::new();
    material.set_density("g/cm3", 2.0).unwrap();
    material.add_nuclide("Li6", 1.0).unwrap();
    let mut nuclide_json_map = std::collections::HashMap::new();
    nuclide_json_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_json_map).unwrap();
    material.calculate_macroscopic_xs(&vec![1], true);
    let material_arc = Arc::new(material);

    let cell1 = Cell::new(Some(1), region1, Some("inner".to_string()), None);
    let cell2 = Cell::new(Some(2), region2.clone(), Some("outer".to_string()), Some(material_arc));

    let geometry = Geometry {
        cells: vec![cell1, cell2],
    };

    // Create source
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 14.06e6,
    };

    let settings = Settings {
        particles: 10000,
        batches: 10,
        source,
        seed: Some(1),
    };

    // Create tallies: one for flux, one for tritium production
    let cell_for_filter = Cell::new(Some(2), region2, Some("outer".to_string()), None);
    let cell_filter = Filter::Cell(CellFilter::new(&cell_for_filter));
    
    let mut flux_tally = Tally::new();
    flux_tally.filters = vec![cell_filter.clone()];
    flux_tally.scores = vec![Score::Flux]; // Use Flux score
    flux_tally.name = Some("flux".to_string());
    flux_tally.initialize_batches(settings.batches as usize);

    // Create energy-binned flux tally
    // Energy bins: logarithmically spaced from 0.01 eV to 20 MeV
    let energy_bins: Vec<f64> = (0..50)
        .map(|i| 10_f64.powf((i as f64) / 49.0 * (20e6_f64.log10() - 0.01_f64.log10()) + 0.01_f64.log10()))
        .collect();
    let energy_filter = EnergyFilter::new(energy_bins);
    
    let mut flux_energy_tally = Tally::new();
    flux_energy_tally.filters = vec![cell_filter, Filter::Energy(energy_filter)];
    flux_energy_tally.scores = vec![Score::Flux];
    flux_energy_tally.name = Some("flux_energy_binned".to_string());
    flux_energy_tally.initialize_batches(settings.batches as usize);

    let tallies = vec![Arc::new(flux_tally), Arc::new(flux_energy_tally)];

    let mut model = Model::new(geometry, settings, tallies);

    let start = Instant::now();
    model.run();
    let elapsed = start.elapsed();

    println!(
        "Simulation completed in {:.6} seconds.",
        elapsed.as_secs_f64()
    );

    // Print results
    let flux_mean = model.tallies[0].get_mean();
    let flux_energy_mean = model.tallies[1].get_mean();
    
    println!("Total Flux: {:.6e}", flux_mean[0]);
    println!("Energy-binned flux has {} bins", flux_energy_mean.len());
    println!("Sum of energy bins: {:.6e}", flux_energy_mean.iter().sum::<f64>());
    
    // Verify that sum of energy bins equals total flux
    let sum_energy_bins: f64 = flux_energy_mean.iter().sum();
    let relative_diff = (sum_energy_bins - flux_mean[0]).abs() / flux_mean[0];
    println!("Relative difference: {:.6e}", relative_diff);
    
    if relative_diff < 0.01 {
        println!("✓ Energy binning is consistent with total flux");
    } else {
        println!("⚠ Warning: Energy bins don't sum to total flux!");
    }
}
