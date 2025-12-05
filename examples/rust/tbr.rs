use yamc::cell::Cell;
use yamc::filter::Filter;
use yamc::filters::CellFilter;
use yamc::geometry::Geometry;
use yamc::model::Model;
use yamc::region::{HalfspaceType, Region, RegionExpr};
use yamc::settings::Settings;
use yamc::source::IndependentSource;
use yamc::stats::AngularDistribution;
use yamc::surface::{BoundaryType, Surface, SurfaceKind};
use yamc::tally::{Score, Tally};
use yamc::Material;
use std::collections::HashMap;
use std::sync::Arc;
use std::time::Instant;

fn main() {
    // Create two-cell geometry: inner sphere (Li6) and outer annular region (Be9)
    let sphere1 = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 1.0,
        },
        boundary_type: BoundaryType::Transmission,
    });

    let sphere2 = Arc::new(Surface {
        surface_id: Some(2),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 200.0,
        },
        boundary_type: BoundaryType::Vacuum,
    });

    // Create regions
    // region1: -sphere1 (inside sphere1)
    let region1 = Region {
        expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
            sphere1.clone(),
        )))),
    };

    // region2: +sphere1 & -sphere2 (between sphere1 and sphere2)
    let region2 = Region {
        expr: RegionExpr::Intersection(
            Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
                sphere1.clone(),
            ))),
            Box::new(RegionExpr::Complement(Box::new(RegionExpr::Halfspace(
                HalfspaceType::Above(sphere2.clone()),
            )))),
        ),
    };

    // Create materials with different absorption characteristics
    let mut material1 = Material::new();
    material1.material_id = Some(1); // Set material_id for MaterialFilter testing

    // Li4SiO4 composition
    material1.add_nuclide("Li6", 0.07 / 2.0).unwrap(); // Li4SiO4
    material1.add_nuclide("Li7", 0.93 / 2.0).unwrap();
    material1.add_nuclide("Be9", 0.5).unwrap();
    // Commented out in original Python:
    // material1.add_element("O", 4.0)
    // material1.add_element("Si", 1.0)
    material1.set_density("g/cm3", 2.0).unwrap();

    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Be9".to_string(), "tests/Be9.json".to_string());
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    nuclide_map.insert("Li7".to_string(), "tests/Li7.json".to_string());
    material1.read_nuclides_from_json(&nuclide_map).unwrap();

    let mat_arc = Arc::new(material1);

    // Create cells
    let cell1 = Cell::new(
        Some(1),
        region1,
        Some("inner_sphere".to_string()),
        None, // No material fill
    );

    let cell2 = Cell::new(
        Some(2),
        region2,
        Some("outer_annular".to_string()),
        Some(mat_arc.clone()),
    );

    let geometry = Geometry::new(vec![cell1.clone(), cell2.clone()]).unwrap();

    // Create source
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 14.06e6,
    };

    let settings = Settings {
        particles: 50000,
        batches: 2,
        source,
        seed: Some(1),
    };

    // Create tallies with CellFilters
    let cell_filter2 = Filter::Cell(CellFilter::new(&cell2));
    let mut tally1 = Tally::new();
    tally1.filters = vec![cell_filter2];
    tally1.scores = vec![Score::MT(105)]; // n,t (tritium production)
    tally1.name = Some("tbr".to_string());

    // Initialize batch data
    tally1.initialize_batches(settings.batches as usize);

    let tallies = vec![Arc::new(tally1)];

    let mut model = Model {
        geometry,
        settings,
        tallies,
    };

    let start = Instant::now();
    model.run();
    let elapsed = start.elapsed();

    println!(
        "Simulation completed in {:.6} seconds.",
        elapsed.as_secs_f64()
    );

    // Tallies are updated in place!
    let mean = model.tallies[0].get_mean();
    println!("TBR (tritium breeding ratio): {}", mean[0]);
}
