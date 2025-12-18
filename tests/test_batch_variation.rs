// Integration test to verify that tally scores show variation between batches
// This ensures the Monte Carlo simulation is properly randomizing and producing
// statistically independent results per batch

use yamc::cell::Cell;
use yamc::tallies::filter::Filter;
use yamc::tallies::CellFilter;
use yamc::geometry::Geometry;
use yamc::model::Model;
use yamc::region::{HalfspaceType, Region, RegionExpr};
use yamc::settings::Settings;
use yamc::source::IndependentSource;
use yamc::distribution_multi::AngularDistribution;
use yamc::surface::{BoundaryType, Surface, SurfaceKind};
use yamc::tallies::tally::{Score, Tally};
use yamc::Material;
use std::collections::HashMap;
use std::sync::Arc;

#[test]
fn test_batch_tally_variation() {
    // Create a sphere as the geometry boundary
    let surf = Arc::new(Surface {
        surface_id: Some(1),
        kind: SurfaceKind::Sphere {
            x0: 0.0,
            y0: 0.0,
            z0: 0.0,
            radius: 200.0,
        },
        boundary_type: BoundaryType::Vacuum,
    });

    // Define region (inside the sphere)
    let region = Region {
        expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(HalfspaceType::Above(
            surf.clone(),
        )))),
    };

    // Create material with Be9 (for MT 16 n,2n reactions)
    let mut material = Material::new();
    material.material_id = Some(1);
    let _ = material.add_nuclide("Be9", 1.0);
    let _ = material.set_density("g/cm3", 1.85);
    let mut nuclide_map = HashMap::new();
    nuclide_map.insert("Be9".to_string(), "tests/Be9.h5".to_string());
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

    // Source: 14 MeV neutrons (high energy for n,2n reactions)
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 14e6,
    };

    // Settings: multiple batches with sufficient particles per batch
    let settings = Settings {
        particles: 1000,
        batches: 10,
        source,
        seed: Some(42),
    };

    // Create tally for MT 16 (n,2n) reactions
    let cell_filter = Filter::Cell(CellFilter::new(&cell));
    let mut tally = Tally::new();
    tally.filters = vec![cell_filter];
    tally.scores = vec![Score::MT(16)]; // MT 16 = (n,2n)
    tally.name = Some("n2n_tally".to_string());

    // Initialize batch data
    tally.initialize_batches(settings.batches as usize);

    let tallies = vec![Arc::new(tally)];

    // Run the model
    let mut model = Model {
        geometry,
        settings,
        tallies,
    };

    model.run();

    // Get the tally results after simulation
    let results = &model.tallies[0];

    // Get batch data for the first (and only) score
    let batch_scores = results.get_batch_data(0);

    println!("\n=== Batch Variation Test Results ===");
    println!("Number of batches: {}", batch_scores.len());
    println!("Batch scores: {:?}", batch_scores);
    println!("Overall mean: {:?}", results.get_mean());
    println!("Overall std dev: {:?}", results.get_std_dev());

    // Check that we have results for all batches
    assert_eq!(
        batch_scores.len(),
        10,
        "Should have tally scores for all 10 batches"
    );

    // Check that there is variation between batches
    // Calculate the range (max - min) and ensure it's non-zero

    let min_score = batch_scores
        .iter()
        .cloned()
        .min_by(|a, b| a.partial_cmp(b).unwrap())
        .expect("Should have at least one batch score");

    let max_score = batch_scores
        .iter()
        .cloned()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .expect("Should have at least one batch score");

    let range = max_score - min_score;
    let mean = results.get_mean()[0];

    println!("Min batch score: {:.6}", min_score);
    println!("Max batch score: {:.6}", max_score);
    println!("Range: {:.6}", range);
    println!("Mean: {:.6}", mean);
    println!("Coefficient of variation: {:.2}%", (range / mean) * 100.0);

    // Assert that there is variation (range > 0)
    assert!(
        range > 0.0,
        "Batch scores should show variation. All batches had identical scores, which suggests \
         the random number generator is not working correctly."
    );

    // Assert that not all batches have the same score
    let all_same = batch_scores.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-10);
    assert!(
        !all_same,
        "Not all batch scores should be identical. This suggests insufficient randomization."
    );

    // Check that the variation is reasonable (not too extreme)
    // For Monte Carlo with ~1000 particles, we expect some variation but not orders of magnitude
    let relative_range = range / mean;
    assert!(
        relative_range > 0.01,
        "Relative variation ({:.4}) is too small - expecting at least 1% variation between batches",
        relative_range
    );

    assert!(
        relative_range < 2.0,
        "Relative variation ({:.4}) is too large - this suggests a problem with the simulation",
        relative_range
    );

    // Check that mean is positive (we should be scoring some reactions)
    assert!(
        mean > 0.0,
        "Mean tally score should be positive - no reactions were tallied"
    );

    // Check that standard deviation is positive (indicates variation)
    let std_dev = results.get_std_dev()[0];
    assert!(
        std_dev > 0.0,
        "Standard deviation should be positive, indicating variation between batches"
    );

    println!("\n✓ Batch variation test passed!");
    println!("  - Batches show variation (range = {:.6})", range);
    println!("  - Relative variation = {:.2}%", relative_range * 100.0);
    println!("  - Standard deviation = {:.6}", std_dev);
}

#[test]
fn test_different_seeds_produce_different_results() {
    // This test verifies that using different seeds produces different results

    fn run_simulation_with_seed(seed: u64) -> Vec<f64> {
        let surf = Arc::new(Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 200.0,
            },
            boundary_type: BoundaryType::Vacuum,
        });

        let region = Region {
            expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(
                HalfspaceType::Above(surf.clone()),
            ))),
        };

        let mut material = Material::new();
        material.material_id = Some(1);
        let _ = material.add_nuclide("Be9", 1.0);
        let _ = material.set_density("g/cm3", 1.85);
        let mut nuclide_map = HashMap::new();
        nuclide_map.insert("Be9".to_string(), "tests/Be9.h5".to_string());
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
            energy: 14e6,
        };

        let settings = Settings {
            particles: 500,
            batches: 5,
            source,
            seed: Some(seed),
        };

        let cell_filter = Filter::Cell(CellFilter::new(&cell));
        let mut tally = Tally::new();
        tally.filters = vec![cell_filter];
        tally.scores = vec![Score::MT(16)];
        tally.name = Some("n2n_tally".to_string());
        tally.initialize_batches(settings.batches as usize);

        let mut model = Model {
            geometry,
            settings,
            tallies: vec![Arc::new(tally)],
        };

        model.run();

        model.tallies[0].get_batch_data(0)
    }

    // Run with two different seeds
    let results_seed1 = run_simulation_with_seed(12345);
    let results_seed2 = run_simulation_with_seed(67890);

    println!("\n=== Different Seeds Test ===");
    println!("Results with seed 12345: {:?}", results_seed1);
    println!("Results with seed 67890: {:?}", results_seed2);

    // Check that the results are different
    let all_different = results_seed1
        .iter()
        .zip(results_seed2.iter())
        .any(|(a, b)| (a - b).abs() > 1e-10);

    assert!(
        all_different,
        "Different seeds should produce different results"
    );

    println!("✓ Different seeds produce different results");
}

#[test]
fn test_same_seed_produces_reproducible_results() {
    // This test verifies that using the same seed produces reproducible results

    fn run_simulation_with_seed(seed: u64) -> Vec<f64> {
        let surf = Arc::new(Surface {
            surface_id: Some(1),
            kind: SurfaceKind::Sphere {
                x0: 0.0,
                y0: 0.0,
                z0: 0.0,
                radius: 200.0,
            },
            boundary_type: BoundaryType::Vacuum,
        });

        let region = Region {
            expr: RegionExpr::Complement(Box::new(RegionExpr::Halfspace(
                HalfspaceType::Above(surf.clone()),
            ))),
        };

        let mut material = Material::new();
        material.material_id = Some(1);
        let _ = material.add_nuclide("Li6", 1.0);
        let _ = material.set_density("g/cm3", 0.5);
        let mut nuclide_map = HashMap::new();
        nuclide_map.insert("Li6".to_string(), "tests/Li6.h5".to_string());
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
            batches: 3,
            source,
            seed: Some(seed),
        };

        let cell_filter = Filter::Cell(CellFilter::new(&cell));
        let mut tally = Tally::new();
        tally.filters = vec![cell_filter];
        tally.scores = vec![Score::MT(105)]; // MT 105 = (n,t)
        tally.name = Some("tritium_tally".to_string());
        tally.initialize_batches(settings.batches as usize);

        let mut model = Model {
            geometry,
            settings,
            tallies: vec![Arc::new(tally)],
        };

        model.run();

        model.tallies[0].get_batch_data(0)
    }

    // Run with the same seed twice
    let results_run1 = run_simulation_with_seed(99999);
    let results_run2 = run_simulation_with_seed(99999);

    println!("\n=== Reproducibility Test ===");
    println!("Results from run 1: {:?}", results_run1);
    println!("Results from run 2: {:?}", results_run2);

    // Check that the results are identical
    assert_eq!(
        results_run1.len(),
        results_run2.len(),
        "Should have same number of batches"
    );

    for (i, (r1, r2)) in results_run1.iter().zip(results_run2.iter()).enumerate() {
        assert!(
            (r1 - r2).abs() < 1e-10,
            "Batch {} results should be identical: {} vs {}",
            i,
            r1,
            r2
        );
    }

    println!("✓ Same seed produces reproducible results");
}
