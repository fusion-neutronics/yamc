use materials_for_mc::*;
use materials_for_mc::model::Model;
use materials_for_mc::source::IndependentSource;
use materials_for_mc::stats::AngularDistribution;
use materials_for_mc::settings::Settings;
use materials_for_mc::tally::Tally;
use materials_for_mc::filter::Filter;
use materials_for_mc::filters::CellFilter;
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
    flux_tally.scores = vec![FLUX_SCORE]; // Use FLUX_SCORE constant
    flux_tally.name = Some("flux".to_string());
    flux_tally.initialize_batches(settings.batches as usize);

    let mut tbr_tally = Tally::new();
    tbr_tally.filters = vec![cell_filter];
    tbr_tally.scores = vec![105]; // n,t (tritium production)
    tbr_tally.name = Some("tbr".to_string());
    tbr_tally.initialize_batches(settings.batches as usize);

    let tallies = vec![Arc::new(flux_tally), Arc::new(tbr_tally)];

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
    let tbr_mean = model.tallies[1].get_mean();
    
    println!("Flux (track-length): {:.6e} cm", flux_mean[0]);
    println!("TBR (tritium breeding ratio): {:.6}", tbr_mean[0]);
}
