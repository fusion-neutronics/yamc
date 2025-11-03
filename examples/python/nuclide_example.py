import materials_for_mc as m4mc
import json

# Configure cross section data paths
m4mc.Config.set_cross_sections({
    "Be9": "tests/Be9.json",
    "Fe54": "tests/Fe54.json",
    "Fe56": "tests/Fe56.json",
    "Fe57": "tests/Fe57.json",
    "Fe58": "tests/Fe58.json",
    "Li6": "tests/Li6.json",
    "Li7": "tests/Li7.json",
    "Pb208": "tests/Pb208.json",  # Added for product examples
})

print("=== Basic Nuclide Properties ===")
nuc1 = m4mc.Nuclide('Li6')
print(f"Nuclide: {nuc1.name}")
print(f"Element: {nuc1.element}")
print(f"Atomic number: {nuc1.atomic_number}")
print(f"Mass number: {nuc1.mass_number}")
print(f"Available temperatures: {nuc1.available_temperatures}")
reaction_mts = nuc1.reaction_mts
if reaction_mts:
    print(f"Available reaction MTs: {reaction_mts[:10]}...")  # Show first 10
else:
    print("No reaction MTs available")

print("\n=== Cross Section Data ===")
# Get cross section data using both MT number and reaction name
micro_n_gamma, energy = nuc1.microscopic_cross_section(reaction='(n,gamma)', temperature='294')
micro_mt_3, energy = nuc1.microscopic_cross_section(reaction=3)

print(f"(n,gamma) cross section at 1 eV: {micro_n_gamma[0]:.6f} barns")
print(f"MT 3 cross section at 1 eV: {micro_mt_3[0]:.6f} barns")
print(f"Energy grid points: {len(energy)}")

print("\n=== Reaction Products and Angular Distributions ===")
# Load a heavier nuclide that has product data
nuc_pb = m4mc.Nuclide('Pb208')
nuc_pb.read_nuclide_from_json('tests/Pb208.json')
reactions = nuc_pb.reactions['294']

# Find reactions with product data
reactions_with_products = []
for mt, reaction in reactions.items():
    reactions_with_products.append((mt, reaction))

print(f"Found {len(reactions_with_products)} reactions with product data")

# Examine the first few reactions with products
for i, (mt, reaction) in enumerate(reactions_with_products[:5]):
    if len(reaction.products) > 0:
        print(f"\n--- Reaction MT {mt} ---")
        product = reaction.products[0]
        print(f"  Particle: {product.particle}")
        print(f"  Emission mode: {product.emission_mode}")
        print(f"  Decay rate: {product.get_decay_rate()} /s")
        print(f"  Number of distributions: {product.num_distributions}")

        # For now, we'll demonstrate sampling instead of examining distributions directly
        # since our Python wrapper doesn't expose the full distribution structure yet
        try:
            # Sample from the reaction product
            outgoing_energy, mu = product.sample(1e6)  # 1 MeV incoming
            print(f"  Sample result at 1 MeV: E_out = {outgoing_energy:.2e} eV, mu = {mu:.3f}")
        except Exception as e:
            print(f"  Sampling failed: {e}")

        continue  # Skip the distribution examination for now
        
        # Note: The following code is commented out because our Python wrapper
        # doesn't currently expose the full distribution structure
        for j, dist in enumerate([]):
            dist_type = type(dist).__name__
            print(f"  Distribution {j+1}: {dist_type}")

            if dist_type == 'UncorrelatedAngleEnergy':
                # This has separate angle and energy distributions
                angle_dist = dist.angle
                print(f"    Angular distribution:")
                print(f"      Energy points: {len(angle_dist.energy)}")
                print(f"      Energy range: {angle_dist.energy[0]:.2e} to {angle_dist.energy[-1]:.2e} eV")
                print(f"      Angle (mu) tabulations: {len(angle_dist.mu)}")

                # Show details of first angular distribution
                if len(angle_dist.mu) > 0:
                    first_mu = angle_dist.mu[0]
                    print(f"      First mu distribution: x = {first_mu.x[:3]}..., p = {first_mu.p[:3]}...")

                # Check if energy distribution is present
                if dist.energy is not None:
                    energy_dist = dist.energy
                    print(f"    Energy distribution: {energy_dist.distribution_type}")

            elif dist_type == 'KalbachMann':
                # This is a correlated angle-energy distribution
                print(f"    Kalbach-Mann distribution:")
                print(f"      Energy points: {len(dist.energy)}")
                print(f"      Energy range: {dist.energy[0]:.2e} to {dist.energy[-1]:.2e} eV")
                print(f"      Energy_out data points: {len(dist.energy_out)}")
                print(f"      Slope data points: {len(dist.slope)}")

            elif dist_type == 'AngleEnergyDistribution':
                # This is a generic wrapper that contains the specific distribution data
                print(f"    Generic AngleEnergyDistribution:")
                print(f"      Distribution type: {dist.distribution_type}")
                # Access the actual data through the data attribute
                # This would contain the specific distribution information

            else:
                print(f"    Unknown distribution type: {dist_type}")
                # Print available attributes for debugging
                attrs = [attr for attr in dir(dist) if not attr.startswith('_')]
                print(f"      Available attributes: {attrs[:5]}...")  # Show first 5 attributes

print("\n=== Reaction Sampling ===")
# Demonstrate reaction sampling (Monte Carlo physics)
print("Sampling reactions at different energies:")
energies_to_test = [1e-3, 1.0, 1e6, 14.1e6]  # thermal, eV, MeV, 14.1 MeV

for energy in energies_to_test:
    reaction = nuc1.sample_reaction(energy=energy, temperature='294', seed=42)
    if reaction:
        print(f"  Energy {energy:.1e} eV -> MT {reaction['mt_number']}")
    else:
        print(f"  Energy {energy:.1e} eV -> No reaction sampled")

print("\n=== Product Sampling Demonstration ===")
# Demonstrate our new sampling functionality with a simple test product
print("Creating and sampling from reaction products:")

# Create test reaction products using our Python API
test_product = m4mc.create_test_reaction_product()
print(f"Test product particle type: {test_product.particle}")

# Sample from it at different energies
for energy in [1e5, 1e6, 5e6]:
    try:
        e_out, mu = test_product.sample(energy)
        print(f"  Input: {energy:.0e} eV -> Output: E={e_out:.2e} eV, mu={mu:.3f}")
    except Exception as e:
        print(f"  Sampling at {energy:.0e} eV failed: {e}")

# Demonstrate sampling functions
print("\nSampling functions:")
mu_isotropic = m4mc.sample_isotropic()
print(f"Isotropic mu sample: {mu_isotropic:.3f}")

# Test tabulated sampling
x_vals = [-1.0, 0.0, 1.0]
p_vals = [0.0, 0.5, 1.0]  # CDF
mu_tabulated = m4mc.sample_tabulated(x_vals, p_vals)
print(f"Tabulated mu sample: {mu_tabulated:.3f}")
