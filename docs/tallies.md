# Tallies

### Scores

Tallies accumulate results during particle transport. Scores can be specified by name or ENDF MT number.

<!-- doctest: skip -->
```python
# Named scores
tally = yamc.Tally(scores=["flux"], cells=cell)
tally = yamc.Tally(scores=["heating"], cells=cell)
tally = yamc.Tally(scores=["heating-local"], cells=cell)     # no secondary photon contribution
tally = yamc.Tally(scores=["absorption"], cells=cell)
tally = yamc.Tally(scores=["fission"], cells=cell)
tally = yamc.Tally(scores=["damage-energy"], cells=cell)     # DPA (MT 444)
tally = yamc.Tally(scores=["H3-production"], cells=cell)     # TBR

# MT numbers
tally = yamc.Tally(scores=[2], cells=cell)    # elastic
tally = yamc.Tally(scores=[18], cells=cell)   # fission
tally = yamc.Tally(scores=[102], cells=cell)  # (n,gamma)
tally = yamc.Tally(scores=[105], cells=cell)  # (n,t) tritium production

# Multiple scores in one tally
tally = yamc.Tally(scores=["flux", "heating", 105], cells=cell)
```

#### Virtual response overlays (``response=``)

By default a tally is *macroscopic*: cross sections are multiplied by the cell's atom density, so it only scores where the responding material physically exists. Pass ``response=`` to instead score a **virtual** response everywhere, even in void -- a detector that isn't present in the geometry. ``response`` and ``nuclides`` are mutually exclusive (``nuclides`` stays the macroscopic per-nuclide breakdown at real in-material density).

The *type* of ``response`` selects the regime:

- ``response="Si28"`` or ``response=["Li6", "Li7"]`` -- microscopic cross sections (barns) at **unit density**; one result bin per nuclide. You weight by atom density yourself in post-processing.
- ``response=<Material>`` -- the material's **macroscopic** response directly, weighting each nuclide's microscopic XS by its real atom density and collapsing to a **single combined bin** (e.g. a SiO2 dose map). No post-processing needed.

<!-- doctest: skip -->
```python
# Nuclide overlay (unit density) -- one bin per nuclide:
tally = yamc.Tally(
    scores=["heating"],
    response="Si28",
    mesh=mesh,
    name="silicon_dose_overlay",
)

# Material overlay -- one combined macroscopic bin, weighted by atom density:
silicon = yamc.Material(composition={"Si28": 1.0}, density=2.33)
tally = yamc.Tally(
    scores=["heating"],
    response=silicon,
    mesh=mesh,
    name="silicon_dose",
)
```

### Estimator

Each tally has an ``estimator`` controlling how a particle history is turned into score contributions:

- ``"track-length"`` (default) -- every cell crossing contributes ``X · weight · Δs``, where ``X`` is the relevant cross section (``Σ_t`` for ``flux`` collapses the XS away into a pure length integral, ``heating_xs`` for ``heating``, ``σ_r`` for a reaction rate, etc.). Lower variance in optically thin regions; touches more bins per history.
- ``"collision"`` -- fires at every collision site with contribution ``(X / Σ_t) · weight``, where ``Σ_t`` is the URR-modified macroscopic total cross section the collision was sampled against. Lower variance in optically thick regions; touches one bin per collision rather than one per crossing. Both estimators converge to the same expected value, so identical seeds must agree to within statistical noise.

<!-- doctest: skip -->
```python
# Default: track-length flux
tally = yamc.Tally(scores=["flux"], cells=cell)

# Collision-estimator flux -- explicit
tally = yamc.Tally(scores=["flux"], cells=cell, estimator="collision")

# Collision-estimator reaction rate (URR-aware for total / elastic /
# fission / absorption / capture; smooth XS otherwise)
tally = yamc.Tally(scores=["(n,gamma)"], cells=cell, estimator="collision")
```

Score / estimator compatibility:

| Score | Track-length | Collision |
| --- | --- | --- |
| ``flux`` | ✓ (default) | ✓ |
| ``heating`` / ``heating-local`` | ✓ -- neutron KERMA (MT 301/901) and photon collision KERMA, both scored along tracks | ✓ -- neutrons use ``(heating_xs / Σ_t) · weight``; photons use the analog absorbed energy ``(E_in − E_out − Σ E_banked_photons) · weight`` at each collision |
| ``reaction-rate`` (``total`` / ``elastic`` / ``fission`` / ``absorption`` / arbitrary MT) | ✓ | ✓ -- URR-aware for MTs 1 / 2 / 18 / 27 / 102 |
| ``H1-production`` / ``H2-production`` / ``H3-production`` / ``He3-production`` / ``He4-production`` | ✓ | ✓ |
| ``damage-energy`` (MT 444) | ✓ | ✓ -- ``(damage_energy_xs / Σ_t) · weight`` per neutron collision |
| ``coherent-scatter`` / ``incoherent-scatter`` / ``photoelectric`` / ``pair-production`` (photon XS) | ✓ | ✓ -- ``(σ_component / Σ_t) · weight`` per photon collision, where ``Σ_t`` is the macroscopic photon total |

Every score type supports both estimators, with these caveats (each refused or documented explicitly rather than silently substituted):

- **Photon heating follows the estimator, with a physical convention split.** ``estimator="collision"`` scores the analog absorbed energy at each photon collision: ``E_in − E_out − Σ E_banked_photons``, electron/positron kinetic energy deposited locally minus the banked TTB bremsstrahlung (the coupled neutron–photon local-deposition convention; matches a reference code's *unfiltered* heating to 0.01% on an absorbing Fe sphere, but not its photon-*filtered* share, which that code books on its electron pseudo-particles). ``estimator="track-length"`` scores collision KERMA along tracks (lower variance, and what mesh tallies want): the same electron energy, but the TTB-radiated part is counted *locally* instead of where those photons land. On thick absorbing bodies the two agree within ~1%; in thin or streaming geometries the track-length value exceeds the absorbed energy by the escaping radiative fraction -- the textbook kerma-vs-dose distinction (4.7% on a 5 MeV broomstick). Photon ``heating-local`` reports the same values as ``heating`` (with no electron transport there is nothing to separate).
- **GPU runs are track-length only.** ``estimator="collision"`` tallies are refused at dispatch with a clear error (issue #353).
- **Unstructured-mesh tallies support both estimators in every mode.** The collision estimator resolves the containing tetrahedron with a BVH point query; under ``woodcock`` / ``hybrid``, flux-score track-length unstructured tallies score true track length along each flight segment (like regular meshes) and other scores use the delta-collision estimator with per-tet attribution.

After running, the chosen estimator is reported on the tally:

<!-- doctest: skip -->
```python
print(tally.estimator)  # "track-length" or "collision"
```

### Variance estimator

Each tally reports a per-bin standard deviation alongside its mean. Variance is estimated end-to-end with per-history Welford accumulation: every source-particle history is one variance sample, per-rayon-worker private scratch is updated online via Welford, then a pairwise reduce folds the workers together. No user-facing knobs -- no batch size, no method selection.

### Spatial binning

Each tally scores in one spatial domain family: cells, materials, mesh, or unstructured mesh. ``cells`` and ``materials`` are mutually exclusive -- a single tally cannot bin over both cells and materials.

<!-- doctest: skip -->
```python
# Cell binning -- single cell (one result bin)
tally = yamc.Tally(scores=["flux"], cells=cell)

# Cell binning -- list of cells (one result bin per cell)
tally = yamc.Tally(scores=["flux"], cells=[cell1, cell2, cell3])
# After results = model.simulate_transport(total_particles=10000):
#   results[tally].mean[0]  → flux in cell1
#   results[tally].mean[1]  → flux in cell2
#   results[tally].mean[2]  → flux in cell3
# A bare single cell is also accepted: cells=cell1 → one bin.

# Material binning -- single material
tally = yamc.Tally(scores=["flux"], materials=material)

# Material binning -- list of materials (one result bin per material)
tally = yamc.Tally(scores=["flux"], materials=[m1, m2])

# Regular mesh binning. `shape` is either per-axis voxel counts [nx, ny, nz]
# or a single int (total voxel count, split into roughly cubic bins over the box).
mesh = yamc.RegularRectangularMesh(
    lower_left=[-10, -10, -10],
    upper_right=[10, 10, 10],
    shape=[20, 20, 20],           # or shape=8000 for ~8000 voxels
)
tally = yamc.Tally(scores=["flux"], mesh=mesh)

# Mesh from geometry bounding box (same shape rules: int total or [nx, ny, nz])
mesh = yamc.RegularRectangularMesh.from_domain(geometry, shape=1_000_000)
tally = yamc.Tally(scores=["flux"], mesh=mesh)

# Regular cylindrical (r, phi, z) mesh -- natural for tokamak / cylindrical
# geometries. phi_bounds defaults to the full circle (0, 2*pi); pass a
# sub-interval (e.g. (0, math.pi)) for a sector.
import math
mesh = yamc.RegularCylindricalMesh(
    r_bounds=(0.0, 10.0),
    z_bounds=(-5.0, 5.0),
    shape=(10, 16, 20),                 # (nr, nphi, nz)
    phi_bounds=(0.0, 2 * math.pi),
    origin=(0.0, 0.0, 0.0),
)
tally = yamc.Tally(scores=["flux"], mesh=mesh)
# Cylindrical cells are NOT equal volume -- outer radial rings are larger.
# Normalise per cell when comparing flux across radii:
#   results[tally].mean[b] / mesh.element_volume(b)

# Unstructured mesh binning (from MeshGeometry)
tally = yamc.Tally(
    scores=["flux"],
    unstructured_mesh=(mesh_geom, "mat:fuel"),  # (MeshGeometry, volume name or ID)
)
```

The ``mean``/``std_dev``/``relative_error`` arrays are laid out flat. When multiple binning axes are present the bin order (slowest → fastest varying) is ``score → cell → material → nuclide → parent_nuclide → energy → mesh``.

### Energy binning

<!-- doctest: skip -->
```python
# Custom bin boundaries (eV)
tally = yamc.Tally(scores=["flux"], cells=cell, energy_bins=[0, 0.625, 1e5, 20e6])

# Named group structure
tally = yamc.Tally(scores=["flux"], cells=cell, energy_group_structure="VITAMIN-J-175")
```

Available group structures: ``"VITAMIN-J-175"`` (175 neutron groups), ``"VITAMIN-J-42"`` (42 photon groups), ``"CCFE-24-PHOTON"`` (24 photon groups).

### Nuclide binning

Score per-nuclide contributions separately.

<!-- doctest: skip -->
```python
tally = yamc.Tally(
    scores=["H3-production"],
    cells=cell,
    nuclides=["Li6", "Li7", "Be9"],
)
# results[tally].mean returns [Li6_score, Li7_score, Be9_score]
```

### Particle binning

Restrict scoring to neutrons or photons. Required for separate neutron/photon spectra.

<!-- doctest: skip -->
```python
tally_n = yamc.Tally(scores=["flux"], cells=cell, particle="neutron")
tally_p = yamc.Tally(scores=["flux"], cells=cell, particle="photon")
```

### Energy function weighting

Apply an energy-dependent weighting function to the tally (e.g. response functions, DPA cross sections).

The values are interpolated with a natural cubic spline, which needs at least
four points. A particle whose energy falls outside the tabulated range scores
nothing into that tally, so make the grid span the whole spectrum you expect.

<!-- doctest: skip -->
```python
tally = yamc.Tally(
    scores=["flux"],
    cells=cell,
    # (energies, values); at least 4 points, energies strictly increasing
    energy_function=([1e-5, 1e3, 1e6, 20e6], [1.0, 1.5, 2.0, 0.5]),
)
```

### Dose coefficients

Apply ICRP-116 dose conversion coefficients directly.

<!-- doctest: skip -->
```python
tally = yamc.Tally(
    scores=["flux"],
    mesh=mesh,
    dose_coefficients=("neutron", "AP"),  # (particle, geometry)
)
```

### Neutron and photon flux spectra

A complete example scoring energy-resolved neutron and photon flux with coupled transport.

<!-- doctest: skip -->
```python
import yamc

yamc.cross_section_data = "endf-b8.1"

# Geometry: iron sphere
sphere = yamc.Sphere(radius=30.0, boundary="vacuum")
iron = yamc.Material(composition={"Fe": 1.0}, density=7.874, temperature=294)
cell = yamc.Cell(region=sphere.below, material=iron)
geometry = yamc.Geometry([cell])

# Source: 14 MeV neutrons
source = yamc.NeutronSource(
    position=(0, 0, 0),
    energy=14.06e6,
)

# Tallies: energy-binned neutron and photon flux
tally_n = yamc.Tally(
    scores=["flux"],
    cells=cell,
    energy_group_structure="VITAMIN-J-175",
    particle="neutron",
    name="neutron_spectrum",
)
tally_p = yamc.Tally(
    scores=["flux"],
    cells=cell,
    energy_group_structure="CCFE-24-PHOTON",
    particle="photon",
    name="photon_spectrum",
)

model = yamc.Model(
    geometry=geometry,
    tallies=[tally_n, tally_p],
    source=source,
    transport_secondary_photons=True,
)
results = model.simulate_transport(total_particles=500000)

# Plot spectra
import numpy as np
import matplotlib.pyplot as plt

n_bins = np.array(tally_n.energy_bins)
n_flux = np.array(results[tally_n].mean)
n_std = np.array(results[tally_n].standard_deviation)

p_bins = np.array(tally_p.energy_bins)
p_flux = np.array(results[tally_p].mean)
p_std = np.array(results[tally_p].standard_deviation)

fig, ax = plt.subplots()
for bins, flux, std, label in [
    (n_bins, n_flux, n_std, "neutron"),
    (p_bins, p_flux, p_std, "photon"),
]:
    centers = np.sqrt(bins[:-1] * bins[1:])      # geometric bin centres
    ax.step(bins[:-1], flux, where="post", label=label)
    ax.errorbar(centers, flux, yerr=std, fmt="none", capsize=2, alpha=0.5)

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("Energy [eV]")
ax.set_ylabel("Flux [cm / source particle]")
ax.legend()
ax.grid(True, alpha=0.3)
fig.savefig("flux_spectra.png", dpi=150)
```


