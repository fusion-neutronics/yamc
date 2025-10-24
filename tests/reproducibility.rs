// Integration test for reproducibility - verifies that simulations with the same seed produce identical results

use std::collections::HashMap;
use std::sync::Arc;
use materials_for_mc::surface::{Surface, SurfaceKind, BoundaryType};
use materials_for_mc::region::{Region, RegionExpr, HalfspaceType};
use materials_for_mc::cell::Cell;
use materials_for_mc::geometry::Geometry;
use materials_for_mc::Material;
use materials_for_mc::tally::Tally;
use materials_for_mc::settings::Settings;
use materials_for_mc::source::IndependentSource;
use materials_for_mc::model::Model;
use materials_for_mc::stats::AngularDistribution;

#[test]
fn test_reproducibility_with_same_seed() {
    // Create geometry
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 2.0 },
        boundary_type: BoundaryType::Vacuum,
    });
    let region = Region { expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(surf.clone())))) };

    // Create material
    let mut material = Material::new();
    material.material_id = Some(1);
    material.add_nuclide("Li6", 1.0);
    material.set_density("g/cm3", 10.0);
    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_map).unwrap();
    material.calculate_macroscopic_xs(&vec![1], true);

    let mat_arc = Arc::new(material);
    let cell = Cell::new(Some(1), region, Some("test_cell".to_string()), Some(mat_arc.clone()));
    let geometry = Geometry::new(vec![cell.clone()]).unwrap();

    // Source
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 1e6,
    };

    // Settings with seed
    let settings = Settings {
        particles: 100,
        batches: 10,
        source: source.clone(),
        seed: Some(42), // Fixed seed for reproducibility
    };

    // Create tally
    let mut tally = Tally::new();
    tally.score = 101; // absorption
    tally.name = Some("test_absorption".to_string());

    // Run simulation 1
    let model1 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![tally.clone()],
    };
    let results1 = model1.run();

    // Run simulation 2 with same seed
    let model2 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![tally.clone()],
    };
    let results2 = model2.run();

    // Run simulation 3 with same seed
    let model3 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![tally.clone()],
    };
    let results3 = model3.run();

    // Verify all three runs produced identical results
    assert_eq!(results1.len(), results2.len(), "Should have same number of tallies");
    assert_eq!(results1.len(), results3.len(), "Should have same number of tallies");

    // Check leakage tally (index 0)
    assert_eq!(results1[0].mean, results2[0].mean, "Leakage should be identical with same seed");
    assert_eq!(results1[0].mean, results3[0].mean, "Leakage should be identical with same seed");
    assert_eq!(results1[0].batch_data, results2[0].batch_data, "Leakage batch data should be identical");
    assert_eq!(results1[0].batch_data, results3[0].batch_data, "Leakage batch data should be identical");

    // Check absorption tally (index 1)
    assert_eq!(results1[1].mean, results2[1].mean, "Absorption should be identical with same seed");
    assert_eq!(results1[1].mean, results3[1].mean, "Absorption should be identical with same seed");
    assert_eq!(results1[1].batch_data, results2[1].batch_data, "Absorption batch data should be identical");
    assert_eq!(results1[1].batch_data, results3[1].batch_data, "Absorption batch data should be identical");

    println!("✓ Reproducibility test passed!");
    println!("  Run 1 - Leakage: {:.6}, Absorption: {:.6}", results1[0].mean, results1[1].mean);
    println!("  Run 2 - Leakage: {:.6}, Absorption: {:.6}", results2[0].mean, results2[1].mean);
    println!("  Run 3 - Leakage: {:.6}, Absorption: {:.6}", results3[0].mean, results3[1].mean);
}

#[test]
fn test_different_seeds_produce_different_results() {
    // Create geometry
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere { x0: 0.0, y0: 0.0, z0: 0.0, radius: 2.0 },
        boundary_type: BoundaryType::Vacuum,
    });
    let region = Region { expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(surf.clone())))) };

    // Create material
    let mut material = Material::new();
    material.material_id = Some(1);
    material.add_nuclide("Li6", 1.0);
    material.set_density("g/cm3", 10.0);
    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_map).unwrap();
    material.calculate_macroscopic_xs(&vec![1], true);

    let mat_arc = Arc::new(material);
    let cell = Cell::new(Some(1), region, Some("test_cell".to_string()), Some(mat_arc.clone()));
    let geometry = Geometry::new(vec![cell.clone()]).unwrap();

    // Source
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 1e6,
    };

    // Settings with different seeds
    let settings1 = Settings {
        particles: 100,
        batches: 10,
        source: source.clone(),
        seed: Some(42),
    };

    let settings2 = Settings {
        particles: 100,
        batches: 10,
        source: source.clone(),
        seed: Some(123),
    };

    // Create tally
    let mut tally = Tally::new();
    tally.score = 101; // absorption
    tally.name = Some("test_absorption".to_string());

    // Run simulation with seed 42
    let model1 = Model {
        geometry: geometry.clone(),
        settings: settings1,
        tallies: vec![tally.clone()],
    };
    let results1 = model1.run();

    // Run simulation with seed 123
    let model2 = Model {
        geometry: geometry.clone(),
        settings: settings2,
        tallies: vec![tally.clone()],
    };
    let results2 = model2.run();

    // Verify different seeds produce different results (with high probability)
    // Note: In principle they could be equal by chance, but with 100 particles this is extremely unlikely
    let different_leakage = results1[0].mean != results2[0].mean;
    let different_absorption = results1[1].mean != results2[1].mean;
    let different_batch_data = results1[1].batch_data != results2[1].batch_data;

    assert!(
        different_leakage || different_absorption || different_batch_data,
        "Different seeds should produce different results (leakage: {} vs {}, absorption: {} vs {})",
        results1[0].mean, results2[0].mean, results1[1].mean, results2[1].mean
    );

    println!("✓ Different seeds test passed!");
    println!("  Seed 42  - Leakage: {:.6}, Absorption: {:.6}", results1[0].mean, results1[1].mean);
    println!("  Seed 123 - Leakage: {:.6}, Absorption: {:.6}", results2[0].mean, results2[1].mean);
}
