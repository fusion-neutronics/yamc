// Integration test for reproducibility - verifies that simulations with the same seed produce identical results

use std::collections::HashMap;
use std::sync::{Arc, Mutex};
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

    let mat_arc = Arc::new(Mutex::new(material));
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

    // Create tallies
    let mut tally1 = Tally::new();
    tally1.score = 101;
    tally1.name = Some("test_absorption_1".to_string());

    let mut tally2 = Tally::new();
    tally2.score = 101;
    tally2.name = Some("test_absorption_2".to_string());

    let mut tally3 = Tally::new();
    tally3.score = 101;
    tally3.name = Some("test_absorption_3".to_string());

    // Run simulation 1
    let mut model1 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![Arc::new(tally1)],
    };
    model1.run();

    // Run simulation 2 with same seed
    let mut model2 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![Arc::new(tally2)],
    };
    model2.run();

    // Run simulation 3 with same seed
    let mut model3 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![Arc::new(tally3)],
    };
    model3.run();

    // Verify all three runs produced identical results
    assert_eq!(model1.tallies.len(), model2.tallies.len(), "Should have same number of tallies");
    assert_eq!(model1.tallies.len(), model3.tallies.len(), "Should have same number of tallies");

    // Check absorption tally (index 0) - no leakage tally anymore
    assert_eq!(model1.tallies[0].get_mean(), model2.tallies[0].get_mean(), "Absorption should be identical with same seed");
    assert_eq!(model1.tallies[0].get_mean(), model3.tallies[0].get_mean(), "Absorption should be identical with same seed");

    // Compare batch data
    let batch_data1 = model1.tallies[0].batch_data.lock().unwrap();
    let batch_data2 = model2.tallies[0].batch_data.lock().unwrap();
    let batch_data3 = model3.tallies[0].batch_data.lock().unwrap();

    for i in 0..batch_data1.len() {
        assert_eq!(batch_data1[i].load(std::sync::atomic::Ordering::SeqCst),
                   batch_data2[i].load(std::sync::atomic::Ordering::SeqCst),
                   "Absorption batch data should be identical");
        assert_eq!(batch_data1[i].load(std::sync::atomic::Ordering::SeqCst),
                   batch_data3[i].load(std::sync::atomic::Ordering::SeqCst),
                   "Absorption batch data should be identical");
    }

    println!("✓ Reproducibility test passed!");
    println!("  Run 1 - Absorption: {:.6}", model1.tallies[0].get_mean());
    println!("  Run 2 - Absorption: {:.6}", model2.tallies[0].get_mean());
    println!("  Run 3 - Absorption: {:.6}", model3.tallies[0].get_mean());
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

    let mat_arc = Arc::new(Mutex::new(material));
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

    // Create tallies
    let mut tally1 = Tally::new();
    tally1.score = 101; // absorption
    tally1.name = Some("test_absorption_1".to_string());

    let mut tally2 = Tally::new();
    tally2.score = 101; // absorption
    tally2.name = Some("test_absorption_2".to_string());

    // Run simulation with seed 42
    let mut model1 = Model {
        geometry: geometry.clone(),
        settings: settings1,
        tallies: vec![Arc::new(tally1)],
    };
    model1.run();

    // Run simulation with seed 123
    let mut model2 = Model {
        geometry: geometry.clone(),
        settings: settings2,
        tallies: vec![Arc::new(tally2)],
    };
    model2.run();

    // Verify different seeds produce different results (with high probability)
    // Note: In principle they could be equal by chance, but with 100 particles this is extremely unlikely
    let different_absorption = model1.tallies[0].get_mean() != model2.tallies[0].get_mean();

    // Check batch data
    let batch_data1 = model1.tallies[0].batch_data.lock().unwrap();
    let batch_data2 = model2.tallies[0].batch_data.lock().unwrap();
    let mut different_batch_data = false;
    println!("Batch data 1: {:?}", batch_data1.iter().map(|a| a.load(std::sync::atomic::Ordering::SeqCst)).collect::<Vec<_>>());
    println!("Batch data 2: {:?}", batch_data2.iter().map(|a| a.load(std::sync::atomic::Ordering::SeqCst)).collect::<Vec<_>>());
    for i in 0..batch_data1.len() {
        if batch_data1[i].load(std::sync::atomic::Ordering::SeqCst) != batch_data2[i].load(std::sync::atomic::Ordering::SeqCst) {
            different_batch_data = true;
            break;
        }
    }

    // For this test, it's possible (though unlikely) that different seeds produce the same mean
    // We'll accept the test passing if either the means OR the batch data are different
    // However, if both are identical, that indicates a problem with seeding
    println!("Different absorption: {}, Different batch data: {}", different_absorption, different_batch_data);

    assert!(
        different_absorption || different_batch_data,
        "Different seeds should produce different results (absorption: {} vs {})",
        model1.tallies[0].get_mean(), model2.tallies[0].get_mean()
    );

    println!("✓ Different seeds test passed!");
    println!("  Seed 42  - Absorption: {:.6}", model1.tallies[0].get_mean());
    println!("  Seed 123 - Absorption: {:.6}", model2.tallies[0].get_mean());
}
