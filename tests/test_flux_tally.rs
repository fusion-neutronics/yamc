// Test flux tallying functionality
use materials_for_mc::cell::Cell;
use materials_for_mc::geometry::Geometry;
use materials_for_mc::model::Model;
use materials_for_mc::region::{HalfspaceType, Region, RegionExpr};
use materials_for_mc::settings::Settings;
use materials_for_mc::source::IndependentSource;
use materials_for_mc::stats::AngularDistribution;
use materials_for_mc::surface::{BoundaryType, Surface, SurfaceKind};
use materials_for_mc::tally::{Tally, Score};
use materials_for_mc::Material;
use std::collections::HashMap;
use std::sync::Arc;

#[test]

#[test]
fn test_score_enum() {
    // Test Score enum conversions
    let flux_score = Score::Flux;
    
    let mt_score = Score::MT(101);
    assert_eq!(mt_score.to_i32(), 101);
    
    // Test from_str
    let parsed_flux = Score::from_str("flux").unwrap();
    
    let invalid_result = Score::from_str("invalid");
    assert!(invalid_result.is_err());
}

#[test]
fn test_flux_tally_simulation() {
    // Create a simple sphere geometry
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 5.0,
        },
        boundary_type: BoundaryType::Vacuum,
    });

    let region = Region {
        expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
            surf.clone(),
        )))),
    };

    // Create material with Li6
    let mut material = Material::new();
    material.material_id = Some(1);
    material.add_nuclide("Li6", 1.0).unwrap();
    material.set_density("g/cm3", 0.46).unwrap();
    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_map).unwrap();

    let mat_arc = Arc::new(material);

    // Create cell
    let cell = Cell::new(
        Some(1),
        region,
        Some("test_cell".to_string()),
        Some(mat_arc.clone()),
    );
    let geometry = Geometry::new(vec![cell.clone()]).unwrap();

    // Create source: point source at origin with 14.1 MeV
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 14.1e6,
    };

    // Create settings
    let settings = Settings {
        particles: 1000,
        batches: 10,
        source,
        seed: Some(42),
    };

    // Create flux tally using Score enum
    let mut flux_tally = Tally::new();
    flux_tally.set_scores_mixed(vec![Score::Flux]);
    flux_tally.name = Some("Flux Tally".to_string());

    // Create model and run
    let mut model = Model {
        geometry,
        settings,
        tallies: vec![Arc::new(flux_tally)],
    };
    model.run();

    // Check results
    let tally_result = &model.tallies[0];
    
    // Flux should be positive
    let mean_flux = tally_result.get_mean()[0];
    assert!(mean_flux > 0.0, "Flux should be positive, got {}", mean_flux);
    
    // Should have proper batch count
    use std::sync::atomic::Ordering;
    assert_eq!(tally_result.n_batches.load(Ordering::Relaxed), 10);
    assert_eq!(tally_result.particles_per_batch.load(Ordering::Relaxed), 1000);
}

#[test]
fn test_flux_and_reaction_mixed_tally() {
    // Create a simple sphere geometry
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 5.0,
        },
        boundary_type: BoundaryType::Vacuum,
    });

    let region = Region {
        expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
            surf.clone(),
        )))),
    };

    // Create material with Li6
    let mut material = Material::new();
    material.material_id = Some(1);
    material.add_nuclide("Li6", 1.0).unwrap();
    material.set_density("g/cm3", 0.46).unwrap();
    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_map).unwrap();

    let mat_arc = Arc::new(material);

    // Create cell
    let cell = Cell::new(
        Some(1),
        region,
        Some("test_cell".to_string()),
        Some(mat_arc.clone()),
    );
    let geometry = Geometry::new(vec![cell.clone()]).unwrap();

    // Create source
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 14.1e6,
    };

    let settings = Settings {
        particles: 1000,
        batches: 10,
        source,
        seed: Some(42),
    };

    // Create tally with both flux and absorption
    let mut tally = Tally::new();
    tally.set_scores_mixed(vec![Score::Flux, Score::MT(101)]);
    tally.name = Some("Mixed Tally".to_string());

    // Create model and run
    let mut model = Model {
        geometry,
        settings,
        tallies: vec![Arc::new(tally)],
    };
    model.run();

    // Check results
    let tally_result = &model.tallies[0];
    
    let means = tally_result.get_mean();
    assert_eq!(means.len(), 2);
    
    // Flux (index 0) should be positive
    assert!(means[0] > 0.0, "Flux should be positive, got {}", means[0]);
    
    // Absorption (index 1) should be non-negative
    assert!(means[1] >= 0.0, "Absorption should be non-negative, got {}", means[1]);
}

#[test]
fn test_set_scores_mixed_reinitializes_storage() {
    // Create tally and set scores
    let mut tally = Tally::new();
    
    // Initial set
    tally.set_scores_mixed(vec![Score::MT(101)]);
    assert_eq!(tally.scores.len(), 1);
    assert_eq!(tally.batch_data.len(), 1);
    assert_eq!(tally.mean.len(), 1);
    assert_eq!(tally.std_dev.len(), 1);
    assert_eq!(tally.rel_error.len(), 1);
    
    // Change scores - should reinitialize all storage
    tally.set_scores_mixed(vec![Score::Flux, Score::MT(101), Score::MT(2)]);
    assert_eq!(tally.scores.len(), 3);
    assert_eq!(tally.batch_data.len(), 3);
    assert_eq!(tally.mean.len(), 3);
    assert_eq!(tally.std_dev.len(), 3);
    assert_eq!(tally.rel_error.len(), 3);
    
    // Verify actual score values
    assert_eq!(tally.scores, vec![Score::Flux, Score::MT(101), Score::MT(2)]);
}
