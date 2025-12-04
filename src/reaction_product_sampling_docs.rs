/*!
# Nuclear Reaction Product Sampling

This module provides comprehensive sampling capabilities for nuclear reaction products,
enabling Monte Carlo neutron transport simulations. The implementation follows the
OpenMC approach for nuclear data sampling.

## Features

- **Angular Distribution Sampling**: Sample scattering cosine (mu) from energy-dependent tabulated distributions
- **Energy Distribution Sampling**: Sample outgoing particle energy from various distribution types
- **Correlated Sampling**: Kalbach-Mann correlated angle-energy sampling for realistic physics
- **Multiple Products**: Handle reactions producing multiple particles (fission, (n,2n), etc.)
- **Distribution Types**: Support for UncorrelatedAngleEnergy, KalbachMann, LevelInelastic, and more

## Usage in Monte Carlo Transport

```rust
use yamc::reaction_product::*;
use rand::thread_rng;

// In your transport loop
let incoming_energy = 14e6; // 14 MeV neutron
let mut rng = thread_rng();

match reaction_mt {
    2 => { // Elastic scattering
        if let Some(product) = reaction.products.get("neutron") {
            let (e_out, mu) = product.sample(incoming_energy, &mut rng);
            particle.energy = e_out;
            particle.scatter_direction(mu, &mut rng);
        }
    },
    16 => { // (n,2n) reaction
        for product in &reaction.products {
            if product.is_particle_type(&ParticleType::Neutron) {
                let (e_out, mu) = product.sample(incoming_energy, &mut rng);
                create_secondary_neutron(e_out, mu);
            }
        }
    },
    18 => { // Fission
        for product in &reaction.products {
            match product.particle {
                ParticleType::Neutron => {
                    // Sample fission neutron
                    let (e_out, mu) = product.sample(incoming_energy, &mut rng);
                    create_fission_neutron(e_out, mu);
                },
                ParticleType::Photon => {
                    // Handle fission gamma rays
                    let (e_out, mu) = product.sample(incoming_energy, &mut rng);
                    create_photon(e_out, mu);
                }
            }
        }
    },
    // Handle other reaction types...
}
```

## Sampling Methods

### 1. ReactionProduct::sample()
Main entry point for sampling outgoing particles:
- Handles multiple distributions with applicability probabilities
- Returns (outgoing_energy, mu_cosine)
- Automatically selects appropriate distribution

### 2. AngleDistribution::sample()
Samples scattering cosine from energy-dependent tabulated distributions:
- Energy interpolation between tabulated points
- CDF inversion for probability sampling
- Fallback to isotropic scattering

### 3. EnergyDistribution::sample()
Samples outgoing energy:
- LevelInelastic: Deterministic energy loss
- Tabulated: Energy-dependent probability distributions
- Extensible to other distribution types

### 4. Kalbach-Mann Sampling
Correlated angle-energy sampling for realistic nuclear physics:
- Samples outgoing energy from tabulated distribution
- Calculates angular distribution parameter from systematics
- Uses exponential distribution f(μ) ∝ exp(a·μ)

## Distribution Types

### UncorrelatedAngleEnergy
Independent sampling of angle and energy:
```rust
let distribution = AngleEnergyDistribution::UncorrelatedAngleEnergy {
    angle: angular_distribution,
    energy: Some(energy_distribution),
};
let (e_out, mu) = distribution.sample(incoming_energy, &mut rng);
```

### KalbachMann
Correlated angle-energy sampling for nuclear reactions:
```rust
let distribution = AngleEnergyDistribution::KalbachMann {
    energy: energy_grid,
    energy_out: outgoing_energy_distributions,
    slope: kalbach_mann_parameters,
};
let (e_out, mu) = distribution.sample(incoming_energy, &mut rng);
```

## Testing

Run the sampling tests to verify implementation:

```bash
# Test basic sampling functionality
cargo run --example test_sampling

# Test with real nuclear data
cargo run --example test_real_sampling
```

## Integration Notes

1. **Energy Units**: All energies are in eV
2. **Angular Units**: mu = cos(θ) where θ is scattering angle
3. **Random Number Generation**: Use thread-safe RNG for parallel transport
4. **Error Handling**: Graceful fallbacks for missing or invalid data
5. **Performance**: Optimized binary search and interpolation routines

## References

- OpenMC Documentation: Nuclear Data Interface
- ENDF-6 Format Manual: Angle-Energy Distributions
- Kalbach-Mann Systematics for Nuclear Reaction Angular Distributions

*/