// Integration test: Elastic scatter only
// This test creates a material with only elastic scattering, no absorption or inelastic
// The source is monoenergetic at 15 MeV, and a flux tally with energy filter is used.
// The expected result is that all neutrons remain at 15 MeV (no energy loss).

use materials_for_mc::*;
use materials_for_mc::model::Model;
use materials_for_mc::source::IndependentSource;
use materials_for_mc::stats::AngularDistribution;
use materials_for_mc::settings::Settings;
use materials_for_mc::tally::{Score, Tally};
use materials_for_mc::filter::Filter;
use materials_for_mc::filters::EnergyFilter;
use std::sync::Arc;

#[test]
fn test_elastic_scatter_only_flux() {
    // Create a material with only elastic scattering
    let mut mat = Material::new("ElasticOnly");
    // Add a single nuclide with only elastic cross section
    let mut nuclide = Nuclide::new("TestNuclide");
    // Set up a constant elastic cross section (e.g., 1 barn) over all energies
    let energy_grid = vec![1e-5, 1e7];
    let xs = vec![1.0, 1.0]; // barns
    nuclide.set_elastic_cross_section(energy_grid.clone(), xs.clone());
    // No absorption, no inelastic, no fission
    mat.add_nuclide(nuclide, 1.0);
    mat.set_density("g/cm3", 1.0);

    // Geometry: single cell filled with this material
    let region = Region::universe();
    let cell = Cell::new(Some(1), region, Some("cell".to_string()), Some(Arc::new(mat)));
    let geometry = Geometry { cells: vec![cell] };

    // Source: monoenergetic at 15 MeV
    let source = IndependentSource {
        space: [0.0, 0.0, 0.0],
        angle: AngularDistribution::Isotropic,
        energy: 15.0e6,
    };

    let settings = Settings {
        particles: 10000,
        batches: 5,
        source,
        seed: Some(42),
    };

    // Energy filter: bins from 1e-5 to 20 MeV
    let energy_bins = vec![1e-5, 20.0e6];
    let energy_filter = EnergyFilter::new(energy_bins);

    // Flux tally with energy filter
    let mut tally = Tally::new();
    tally.filters = vec![Filter::Energy(energy_filter)];
    tally.scores = vec![Score::Flux];
    tally.name = Some("flux_energy_binned".to_string());
    tally.initialize_batches(settings.batches as usize);

    let tallies = vec![Arc::new(tally)];
    let mut model = Model::new(geometry, settings, tallies);
    model.run();

    // Get the flux tally result
    let flux = model.tallies[0].get_mean();
    // For realistic elastic kinematics, neutrons will lose energy with each scatter (except for n-H)
    // So the flux should be spread below the source energy, not all at 15 MeV
    // With only one wide bin, all flux is in that bin, but let's use more bins to check the spectrum
    println!("Elastic-only flux tally: {:?}", flux);
    // The highest energy bin should have nonzero flux, but lower bins should also have some flux
    let n_bins = flux.len();
    assert!(n_bins > 1, "Test should use multiple energy bins to check spectrum");
    assert!(flux[n_bins - 1] > 0.0, "Highest energy bin should have flux");
    let lower_flux: f64 = flux[..n_bins - 1].iter().sum();
    assert!(lower_flux > 0.0, "Lower energy bins should have some flux due to downscatter");
}
