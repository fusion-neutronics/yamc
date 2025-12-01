import time
import materials_for_mc as mc
import matplotlib.pyplot as plt
import numpy as np
    
# Create two-cell geometry: inner sphere (Li6) and outer annular region (Be9)
sphere1 = mc.Sphere(
    surface_id=1,
    x0=0.0,
    y0=0.0,
    z0=0.0,
    r=1.0,
    boundary_type='transmission',
)
sphere2 = mc.Sphere(
    surface_id=2,
    x0=0.0,
    y0=0.0,
    z0=0.0,
    r=200.0,
    boundary_type='vacuum',
)
region1 = -sphere1
region2 = +sphere1 & -sphere2

# Create materials with different absorption characteristics
material1 = mc.Material()
material1.material_id = 1  # Set material_id for MaterialFilter testing
material1.add_nuclide("Li6", 0.07/2)  # Li4SiO4
material1.add_nuclide("Li7", 0.93/2)
material1.add_nuclide("Be9", 0.5)
# material1.add_element("O", 4.0)
# material1.add_element("Si", 1.0)
material1.set_density("g/cm3", 2.0)
material1.read_nuclides_from_json({"Be9": "tests/Be9.json", "Li6": "tests/Li6.json", "Li7": "tests/Li7.json"})


# Create cells
cell1 = mc.Cell(
    cell_id=1,
    name="inner_sphere",
    region=region1,
)
cell2 = mc.Cell(
    cell_id=2,
    name="outer_annular",
    region=region2,
    fill=material1,
)
geometry = mc.Geometry([cell1, cell2])

source = mc.IndependentSource(
    space=[0.0, 0.0, 0.0],
    angle=mc.stats.Isotropic(),
    energy=1.06e6
)
settings = mc.Settings(
    particles=50000,
    batches=2,
    source=source,
    seed=1
)

# Create tallies with CellFilters
cell_filter2 = mc.CellFilter(cell2)
tally1 = mc.Tally()
tally1.filters = [cell_filter2]
tally1.scores = ['flux']
tally1.name = "flux"

# Create energy-binned flux tally
# Energy bins: logarithmically spaced from 0.1 eV to 20 MeV
energy_bins = np.logspace(np.log10(0.1), np.log10(20e6), 50)  # 49 bins in eV
energy_filter = mc.EnergyFilter(energy_bins)

tally2 = mc.Tally()
tally2.filters = [cell_filter2, energy_filter]
tally2.scores = ['flux']
tally2.name = "flux_energy_binned"

tallies = [tally1, tally2]

model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)

time.start = time.time()
model.run()
print(f"Simulation completed in {time.time() - time.start} seconds.")
# Tallies are updated in place!
print(f"Total Flux: {tally1.mean}")
print(f"Energy-binned flux has {len(tally2.mean)} bins")

# Plot the energy-binned flux spectrum
fig, ax = plt.subplots(figsize=(10, 6))

# Get flux values from energy-binned tally
# tally2.mean has one value per energy bin
flux_values = tally2.mean

# Plot as step function
ax.step(energy_bins[:-1], flux_values, where='post', label='Flux per energy bin')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Energy (eV)')
ax.set_ylabel('Flux (particles/cm²/source particle)')
ax.set_title('Energy-dependent Flux Spectrum')
ax.grid(True, alpha=0.3)
ax.legend()

plt.tight_layout()
plt.savefig('flux_spectrum.png', dpi=150)
print("Flux spectrum plot saved as 'flux_spectrum.png'")
plt.show()
