// Integration test for reproducibility - verifies that simulations with the same seed produce identical results

use materials_for_mc::cell::Cell;
use materials_for_mc::geometry::Geometry;
use materials_for_mc::model::Model;
use materials_for_mc::region::{HalfspaceType, Region, RegionExpr};
use materials_for_mc::settings::Settings;
use materials_for_mc::source::IndependentSource;
use materials_for_mc::stats::AngularDistribution;
use materials_for_mc::surface::{BoundaryType, Surface, SurfaceKind};
use materials_for_mc::tally::{Score, Tally};
use materials_for_mc::Material;
use std::collections::HashMap;
use std::sync::Arc;

#[test]
fn test_separate_vs_multi_score_tallies_equivalence() {
    // Setup: geometry, material, source, settings
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 2.0,
        },
        boundary_type: BoundaryType::Vacuum,
    });
    let region = Region {
        expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
            surf.clone(),
        )))),
    };

    let mut material = Material::new();
    material.material_id = Some(1);
    let _ = material.add_nuclide("Li6", 1.0);
    let _ = material.set_density("g/cm3", 10.0);
    let mut nuclide_map = std::collections::HashMap::new();
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_map).unwrap();
    let mat_arc = Arc::new(material);
    let cell = Cell::new(
        Some(1),
        region,
        Some("test_cell".to_string()),
        Some(mat_arc.clone()),
    );
    let geometry = Geometry::new(vec![cell.clone()]).unwrap();

    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 1e6,
    };
    let settings = Settings {
        particles: 100,
        batches: 10,
        source: source.clone(),
        seed: Some(12345),
    };
    let num_batches = settings.batches as usize;

    // --- Two tallies, one score each ---
    let mut tally_a = Tally::new();
    tally_a.scores = vec![Score::MT(101)];
    tally_a.name = Some("absorption_101".to_string());
    tally_a.initialize_batches(num_batches);

    let mut tally_b = Tally::new();
    tally_b.scores = vec![Score::MT(102)];
    tally_b.name = Some("absorption_102".to_string());
    tally_b.initialize_batches(num_batches);

    let mut model_sep = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![Arc::new(tally_a), Arc::new(tally_b)],
    };
    model_sep.run();
    let mean_a = model_sep.tallies[0].get_mean()[0];
    let mean_b = model_sep.tallies[1].get_mean()[0];

    // --- One tally, two scores ---
    let mut tally_multi = Tally::new();
    tally_multi.scores = vec![Score::MT(101), Score::MT(102)];
    tally_multi.name = Some("absorption_101_102".to_string());
    tally_multi.initialize_batches(num_batches);

    let mut model_multi = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![Arc::new(tally_multi)],
    };
    model_multi.run();
    let means_multi = &model_multi.tallies[0].get_mean();

    assert!(
        (mean_a - means_multi[0]).abs() < 1e-12,
        "Separate tally 101 and multi-score tally[0] should match"
    );
    assert!(
        (mean_b - means_multi[1]).abs() < 1e-12,
        "Separate tally 102 and multi-score tally[1] should match"
    );

    // --- Three tallies, one score each ---

    // Create fresh tallies for three-tally test
    let mut tally_a3 = Tally::new();
    tally_a3.scores = vec![Score::MT(101)];
    tally_a3.name = Some("absorption_101".to_string());
    tally_a3.initialize_batches(num_batches);

    let mut tally_b3 = Tally::new();
    tally_b3.scores = vec![Score::MT(102)];
    tally_b3.name = Some("absorption_102".to_string());
    tally_b3.initialize_batches(num_batches);

    let mut tally_c = Tally::new();
    tally_c.scores = vec![Score::MT(103)];
    tally_c.name = Some("absorption_103".to_string());
    tally_c.initialize_batches(num_batches);

    let mut model_sep3 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![Arc::new(tally_a3), Arc::new(tally_b3), Arc::new(tally_c)],
    };
    model_sep3.run();
    let mean_a3 = model_sep3.tallies[0].get_mean()[0];
    let mean_b3 = model_sep3.tallies[1].get_mean()[0];
    let mean_c3 = model_sep3.tallies[2].get_mean()[0];

    // --- One tally, three scores ---
    let mut tally_multi3 = Tally::new();
    tally_multi3.scores = vec![Score::MT(101), Score::MT(102), Score::MT(103)];
    tally_multi3.name = Some("absorption_101_102_103".to_string());
    tally_multi3.initialize_batches(num_batches);

    let mut model_multi3 = Model {
        geometry: geometry.clone(),
        settings: settings.clone(),
        tallies: vec![Arc::new(tally_multi3)],
    };
    model_multi3.run();
    let means_multi3 = &model_multi3.tallies[0].get_mean();

    assert!(
        (mean_a3 - means_multi3[0]).abs() < 1e-12,
        "Separate tally 101 and multi-score tally[0] should match"
    );
    assert!(
        (mean_b3 - means_multi3[1]).abs() < 1e-12,
        "Separate tally 102 and multi-score tally[1] should match"
    );
    assert!(
        (mean_c3 - means_multi3[2]).abs() < 1e-12,
        "Separate tally 103 and multi-score tally[2] should match"
    );

    println!("✓ Separate vs multi-score tally equivalence test passed!");
}

#[test]
fn test_reproducibility_with_same_seed() {
    // Create geometry
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 2.0,
        },
        boundary_type: BoundaryType::Vacuum,
    });
    let region = Region {
        expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
            surf.clone(),
        )))),
    };

    // Create material
    let mut material = Material::new();
    material.material_id = Some(1);
    let _ = material.add_nuclide("Li6", 1.0);
    let _ = material.set_density("g/cm3", 10.0);
    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_map).unwrap();

    let mat_arc = Arc::new(material);
    let cell = Cell::new(
        Some(1),
        region,
        Some("test_cell".to_string()),
        Some(mat_arc.clone()),
    );
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
    tally1.scores = vec![Score::MT(101)];
    tally1.name = Some("test_absorption_1".to_string());
    let mut tally2 = Tally::new();
    tally2.scores = vec![Score::MT(101)];
    tally2.name = Some("test_absorption_2".to_string());
    let mut tally3 = Tally::new();
    tally3.scores = vec![Score::MT(101)];
    tally3.name = Some("test_absorption_3".to_string());
    // Initialize batch data for all tallies before wrapping in Arc
    let num_batches = settings.batches as usize;
    tally1.initialize_batches(num_batches);
    tally2.initialize_batches(num_batches);
    tally3.initialize_batches(num_batches);

    // Run simulation 1
    let mut model1 = Model::new(
        geometry.clone(),
        settings.clone(),
        vec![Arc::new(tally1)],
    );
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
    assert_eq!(
        model1.tallies.len(),
        model2.tallies.len(),
        "Should have same number of tallies"
    );
    assert_eq!(
        model1.tallies.len(),
        model3.tallies.len(),
        "Should have same number of tallies"
    );

    // Check absorption tally (index 0) - no leakage tally anymore
    assert_eq!(
        model1.tallies[0].get_mean()[0],
        model2.tallies[0].get_mean()[0],
        "Absorption should be identical with same seed"
    );
    assert_eq!(
        model1.tallies[0].get_mean()[0],
        model3.tallies[0].get_mean()[0],
        "Absorption should be identical with same seed"
    );

    // Compare batch data
    use std::sync::atomic::Ordering;
    for i in 0..model1.tallies[0].batch_data.len() {
        let v1 = model1.tallies[0].batch_data[i].lock().unwrap();
        let v2 = model2.tallies[0].batch_data[i].lock().unwrap();
        let v3 = model3.tallies[0].batch_data[i].lock().unwrap();
        let vals1: Vec<u64> = v1.iter().map(|a| a.load(Ordering::SeqCst)).collect();
        let vals2: Vec<u64> = v2.iter().map(|a| a.load(Ordering::SeqCst)).collect();
        let vals3: Vec<u64> = v3.iter().map(|a| a.load(Ordering::SeqCst)).collect();
        assert_eq!(vals1, vals2, "Absorption batch data should be identical");
        assert_eq!(vals1, vals3, "Absorption batch data should be identical");
    }

    println!("✓ Reproducibility test passed!");
    println!(
        "  Run 1 - Absorption: {:?}",
        model1.tallies[0].get_mean()[0]
    );
    println!(
        "  Run 2 - Absorption: {:?}",
        model2.tallies[0].get_mean()[0]
    );
    println!(
        "  Run 3 - Absorption: {:?}",
        model3.tallies[0].get_mean()[0]
    );
}

#[test]
fn test_different_seeds_produce_different_results() {
    // Create geometry
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 2.0,
        },
        boundary_type: BoundaryType::Vacuum,
    });
    let region = Region {
        expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
            surf.clone(),
        )))),
    };

    // Create material
    let mut material = Material::new();
    material.material_id = Some(1);
    let _ = material.add_nuclide("Li6", 1.0);
    let _ = material.set_density("g/cm3", 10.0);
    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Li6".to_string(), "tests/Li6.json".to_string());
    material.read_nuclides_from_json(&nuclide_map).unwrap();

    let mat_arc = Arc::new(material);
    let cell = Cell::new(
        Some(1),
        region,
        Some("test_cell".to_string()),
        Some(mat_arc.clone()),
    );
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
    tally1.scores = vec![Score::MT(101)]; // absorption
    tally1.name = Some("test_absorption_1".to_string());
    let mut tally2 = Tally::new();
    tally2.scores = vec![Score::MT(101)]; // absorption
    tally2.name = Some("test_absorption_2".to_string());
    // Initialize batch data for all tallies before wrapping in Arc
    let num_batches = settings1.batches as usize;
    tally1.initialize_batches(num_batches);
    tally2.initialize_batches(num_batches);

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
    let different_absorption = model1.tallies[0].get_mean()[0] != model2.tallies[0].get_mean()[0];

    // Check batch data
    let mut different_batch_data = false;
    use std::sync::atomic::Ordering;
    for i in 0..model1.tallies[0].batch_data.len() {
        let v1 = model1.tallies[0].batch_data[i].lock().unwrap();
        let v2 = model2.tallies[0].batch_data[i].lock().unwrap();
        let vals1: Vec<u64> = v1.iter().map(|a| a.load(Ordering::SeqCst)).collect();
        let vals2: Vec<u64> = v2.iter().map(|a| a.load(Ordering::SeqCst)).collect();
        println!("Batch data 1[{}]: {:?}", i, vals1);
        println!("Batch data 2[{}]: {:?}", i, vals2);
        if vals1 != vals2 {
            different_batch_data = true;
            break;
        }
    }

    // For this test, it's possible (though unlikely) that different seeds produce the same mean
    // We'll accept the test passing if either the means OR the batch data are different
    // However, if both are identical, that indicates a problem with seeding
    println!(
        "Different absorption: {}, Different batch data: {}",
        different_absorption, different_batch_data
    );

    assert!(
        different_absorption || different_batch_data,
        "Different seeds should produce different results (absorption: {:?} vs {:?})",
        model1.tallies[0].get_mean()[0],
        model2.tallies[0].get_mean()[0]
    );

    println!("✓ Different seeds test passed!");
    println!(
        "  Seed 42  - Absorption: {:?}",
        model1.tallies[0].get_mean()[0]
    );
    println!(
        "  Seed 123 - Absorption: {:?}",
        model2.tallies[0].get_mean()[0]
    );
}
