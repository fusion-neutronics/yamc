# Sources

### Source

A source defines where, in what direction, and at what energy particles are born. All parameters are keyword-only.

```python
source = yamc.NeutronSource(
    position=(0, 0, 0),
    energy=14.06e6,
)
```

Defaults: origin at ``[0, 0, 0]``, 14.06 MeV, isotropic, neutron, strength 1.0.

``position`` accepts a plain ``(x, y, z)`` tuple and ``energy`` a single value in
eV as shorthand for a fixed point and a monoenergetic source. Pass an explicit
spatial or energy distribution (see below) when you need a ring, a spectrum, or
any other multi-location / multi-energy source.

The spatial, energy and angular distributions below live in the ``yamc.sources``
submodule (``yamc.sources.Discrete``, ``yamc.sources.Isotropic``, ...). The
``NeutronSource`` / ``PhotonSource`` builders themselves stay at the top level;
you can also ``from yamc import sources`` and write ``sources.Discrete(...)``.

### Spatial distributions

```python
# Fixed point, using the (x, y, z) tuple shorthand
source = yamc.NeutronSource(position=(0, 0, 0), energy=14.06e6)

# Cylindrical ring source (a bare ring). For a real tokamak plasma, with
# its shape and its density/temperature profiles, see
# "Parametric tokamak plasma source" below.
import math
source = yamc.NeutronSource(
    position=yamc.sources.CylindricalRing(
        radius=yamc.sources.Discrete([300], [1]),   # radius in cm
        phi=yamc.sources.Uniform(0, 2 * math.pi),   # full ring
        z=yamc.sources.Discrete([0], [1]),          # midplane
    ),
    energy=14.06e6,
)
```

### Source on a mesh

``yamc.sources.MeshDistribution`` spreads a source over the elements of a mesh:
the voxels of a ``RegularRectangularMesh``, the cells of a
``RegularCylindricalMesh``, or the tetrahedra of a tet ``MeshGeometry``. Each
element emits in proportion to its strength, from a point uniform within it.
Pass it as ``position=`` to ``NeutronSource`` or ``PhotonSource``.

```python
mesh = yamc.RegularRectangularMesh([0, 0, 0], [10, 10, 10], [2, 1, 1])
space = yamc.sources.MeshDistribution(mesh, strengths=[3.0, 1.0])
source = yamc.PhotonSource(position=space, energy=1e6)

print(space.num_elements)    # 2
print(space.total_strength)  # 4.0
```

Here three quarters of the photons are born in the voxel with ``x`` from 0 to
5 cm and one quarter in the other.

``strengths`` takes one non-negative value per element, in the order a tally
on the same mesh numbers its bins, so a mesh tally's ``mean`` can be passed in
directly. Omitted, each element's strength is its volume, which makes the
source uniform over the mesh. Elements with zero strength never emit. A list of
the wrong length, a NaN, infinite or negative value, or strengths that sum to
zero raise ``ValueError``. When the strengths are rates (particles per second),
``total_strength`` is the source rate to scale per-source-particle tallies by.

``cells=`` (a ``Cell`` or a list of them) or ``materials=`` (a ``Material`` or
a list) keeps each element's strength but has it emit only from the parts of it
inside those cells, or inside the cells filled with those materials. A decay
photon source on a voxel half filled with steel is then born in the steel, not
in the void beside it. These need the geometry, so they apply in a ``Model``
run; ``space.sample()`` on its own draws from the whole element.

### Energy distributions

```python
# Monoenergetic
energy = 14.06e6

# Discrete line spectrum, e.g. Co-60's two main decay gammas,
# 1.17 MeV and 1.33 MeV, emitted in ~equal intensity
energy = yamc.sources.Discrete([1.1732e6, 1.3325e6], [0.5, 0.5])

# Histogram (piecewise-constant) spectrum: n+1 ascending bin edges (eV) and one
# per-bin weight. Weights are per-bin probability mass and are normalized
# internally, so a raw multigroup flux per group can be passed directly. Here
# 30% of neutrons are born uniformly in [0, 1] MeV and 70% in [1, 20] MeV.
energy = yamc.sources.Histogram([0.0, 1e6, 20e6], [0.3, 0.7])

# The bin edges can also be a built-in group structure name instead of a list
# (e.g. "CCFE-709", "UKAEA-1102"); pass one weight per group. yamc.group_structure()
# returns the edges of any of them, and yamc.group_structure_names() lists them;
# see the table under Tallies > Energy binning.
energy = yamc.sources.Histogram("CCFE-709", flux_per_group)   # 709 weights

# Uniform range
energy = yamc.sources.Uniform(1e6, 20e6)

# Gaussian
energy = yamc.sources.Normal(14.06e6, 3e5)

# DT fusion neutron spectrum (Ballabio et al.)
energy = yamc.sources.fusion_neutron_spectrum(ion_temp=20000.0)
# DD fusion neutron spectrum
energy = yamc.sources.fusion_neutron_spectrum(ion_temp=20000.0, reactants='DD')
```

### Angular distributions

```python
direction = yamc.sources.Isotropic()                 # default, usually omitted
direction = yamc.sources.Monodirectional([0, 0, 1])  # beam in +z
```

### Parametric tokamak plasma source

`yamc.sources.tokamak_source()` builds a tokamak's neutron source from the
plasma itself (its shape, confinement mode, and ion density and temperature
profiles) rather than from a hand-placed ring.

```python
import yamc

# Shape (cm, dimensionless)
shape = yamc.sources.TokamakShape(
    major_radius=906.0,
    minor_radius=292.258,
    elongation=1.557,
    triangularity=0.270,
    shafranov_factor=0.44789,
)

plasma = dict(
    shape=shape,
    # Confinement mode: "L", "H" or "A"
    mode="H",
    # Ion density profile (m^-3, cm for the pedestal radius)
    ion_density_centre=1.09e20,
    ion_density_peaking_factor=1,
    ion_density_pedestal=1.09e20,
    ion_density_separatrix=3e19,
    pedestal_radius=0.8 * 292.258,
    # Ion temperature profile (eV)
    ion_temperature_centre=45.9e3,
    ion_temperature_peaking_factor=8.06,
    ion_temperature_beta=6.0,
    ion_temperature_pedestal=6.09e3,
    ion_temperature_separatrix=0.1e3,
)

sources = yamc.sources.tokamak_source(**plasma)
model = yamc.Model(geometry=geometry, source=sources, tallies=tallies)
```

Lengths are in cm, ion densities in m^-3, ion temperatures in eV and angles in
radians. The full argument list is in the [API reference](api.md).

The shape is its own object, `TokamakShape`, so it can be checked before any
source is built. It refuses an unphysical shape (a minor radius not less than
the major radius, a triangularity outside -1 to 1, a Shafranov factor of half
the minor radius or more) when it is constructed, and gives the plasma volume
and points on any flux surface:

```python
shape.volume()                  # plasma volume inside the boundary, cm^3
r, z = shape.evaluate(1.0, 0.0) # outboard midplane of the boundary, cm
```

`evaluate(rho, theta)` takes the normalised radius `rho` (0 on the magnetic
axis, 1 on the boundary) and the poloidal angle `theta` in radians. Either may
be a sequence, so a whole flux surface comes back in one call.

It returns a list of `NeutronSource` ring sources. The plasma is
discretised onto a cylindrical (R, phi, Z) mesh and each non-empty voxel
becomes one source per reaction, carrying the Ballabio spectrum at that voxel's
emission-weighted ion temperature. Strengths are normalised to sum to 1, so the
absolute neutron rate comes from the simulation's source rate. Pass the list
straight to `Model(source=...)`, as under
"Multiple sources and strength weighting" below.

The mesh defaults to `(100, 100)`, which is about 15000 sources for D-T fuel.
Every non-empty voxel costs one source, so pass a coarser
`mesh_resolution=(50, 50)` when a smoother source is not worth the memory.

#### Toroidal sectors

A wedge model only needs the matching wedge of plasma. `start_angle` and
`rotation_angle` are in radians; a negative `rotation_angle` sweeps the other
way.

```python
import math

sources = yamc.sources.tokamak_source(
    **plasma,
    start_angle=0.0,
    rotation_angle=math.pi / 4,   # a 45 degree sector
)
```

#### Fuel

`fuel` takes atom fractions that sum to 1, and defaults to 50:50 D-T. The
supported species are `"D"` and `"T"`, giving D(d,n)3He and T(d,n)4He neutrons.

```python
sources = yamc.sources.tokamak_source(**plasma, fuel={"D": 1.0})          # D-D only
sources = yamc.sources.tokamak_source(**plasma, fuel={"D": 0.4, "T": 0.6})
```

T(t,2n)4He is not included because yamc has no T-T neutron spectrum model, so
D-D and D-T carry the whole yield. For 50:50 D-T fuel at reactor temperatures
the T-T channel is well under a percent of the total (0.7% for the parameters
above), and pure-tritium fuel is rejected rather than silently returning no
neutrons.

#### Plasma profiles

The profiles the source is built from are exposed on their own, for plotting
and for sanity-checking a parameter set. They depend only on the minor radius,
and each accepts a float or a sequence of floats.

```python
densities = yamc.sources.tokamak_ion_density(
    mode="H",
    ion_density_centre=1.09e20,
    ion_density_peaking_factor=1,
    ion_density_pedestal=1.09e20,
    ion_density_separatrix=3e19,
    minor_radius=292.258,
    pedestal_radius=0.8 * 292.258,
    r=[0.0, 100.0, 292.258],           # cm
)

temperatures = yamc.sources.tokamak_ion_temperature(
    r=[0.0, 100.0, 292.258],
    mode="H",
    ion_temperature_centre=45.9e3,
    ion_temperature_peaking_factor=8.06,
    ion_temperature_beta=6.0,
    ion_temperature_pedestal=6.09e3,
    ion_temperature_separatrix=0.1e3,
    minor_radius=292.258,
    pedestal_radius=0.8 * 292.258,
)
```

`tokamak_neutron_source_density()` completes the set: it gives the Bosch-Hale
reaction rate density in neutrons/s/m^3. Mapping flux coordinates onto
`(R, Z)`, with the elongation, triangularity and Shafranov shift applied, is
`TokamakShape.evaluate()` above.

Plotting the two profiles across the plasma:

```python
import matplotlib.pyplot as plt
import yamc

minor_radius = 292.258
profile = dict(mode="H", minor_radius=minor_radius, pedestal_radius=0.8 * minor_radius)
radii = [i * minor_radius / 200 for i in range(201)]

densities = yamc.sources.tokamak_ion_density(
    **profile,
    ion_density_centre=1.09e20,
    ion_density_peaking_factor=1,
    ion_density_pedestal=1.09e20,
    ion_density_separatrix=3e19,
    r=radii,
)
temperatures = yamc.sources.tokamak_ion_temperature(
    **profile,
    r=radii,
    ion_temperature_centre=45.9e3,
    ion_temperature_peaking_factor=8.06,
    ion_temperature_beta=6.0,
    ion_temperature_pedestal=6.09e3,
    ion_temperature_separatrix=0.1e3,
)

fig, (top, bottom) = plt.subplots(2, 1, sharex=True, figsize=(8, 6))
top.plot(radii, densities)
top.set_ylabel("Ion density [m$^{-3}$]")
bottom.plot(radii, temperatures)
bottom.set_xlabel("Minor radius [cm]")
bottom.set_ylabel("Ion temperature [eV]")
for axis in (top, bottom):
    axis.axvline(0.8 * minor_radius, ls="--", c="grey", lw=1, label="pedestal")
    axis.grid(True, alpha=0.3)
top.legend()
plt.show()
```

#### Plotting the emission

Sampling the sources shows the plasma's poloidal cross-section: the D-shape
from the triangularity, the elongation, and the emission concentrated near the
magnetic axis.

```python
import math

import matplotlib.pyplot as plt
import numpy as np
import yamc

sources = yamc.sources.tokamak_source(**plasma, mesh_resolution=(60, 60))

# Share the points out between the sources by strength, the way the transport
# picks them. Sampling every source equally would instead give each voxel the
# same weight, hiding both the peaked emission and the true D-D fraction.
n_points = 40_000
strengths = np.array([source.strength for source in sources])
counts = np.random.default_rng(1).multinomial(n_points, strengths / strengths.sum())

radii, heights, energies = [], [], []
for source, n in zip(sources, counts):
    if n == 0:
        continue
    positions, sampled = source.sample_n(int(n))
    for (x, y, z), energy in zip(positions, sampled):
        radii.append(math.hypot(x, y))
        heights.append(z)
        energies.append(energy / 1e6)

fig, ax = plt.subplots(figsize=(7, 7))
scatter = ax.scatter(radii, heights, c=energies, s=2, alpha=0.4, cmap="plasma")
ax.set_xlabel("R [cm]")
ax.set_ylabel("Z [cm]")
ax.set_aspect("equal")
ax.set_title("Neutron birth sites, coloured by energy")
fig.colorbar(scatter, ax=ax, label="Birth energy [MeV]")
plt.show()
```

Emission crowds towards the magnetic axis, where the plasma is hottest, and
the sparse 2.5 MeV points are the D-D component, about 0.5% of the yield for
50:50 D-T fuel at these temperatures.

### Multiple sources and strength weighting

Pass a list of sources with relative strengths to ``Model()``.

<!-- doctest: skip -->
```python
source_a = yamc.NeutronSource(
    position=(0, 0, -7),
    energy=14.06e6,
    strength=2.0,
)
source_b = yamc.NeutronSource(
    position=(0, 0, 7),
    energy=2.5e6,
    strength=1.0,
)
model = yamc.Model(geometry=geometry, source=[source_a, source_b])
```

### VTK-HDF source export

Sampled source points can be exported for [ParaView](https://www.paraview.org/) visualization.

<!-- doctest: skip -->
```python
source.to_vtkhdf("source_points.vtkhdf", samples=1000)
```

The viewer below is rendered live when these docs are built: it samples a tokamak-like cylindrical-ring plasma source, writes ``source_points.vtkhdf``, reads it back, and shows the sampled emission sites as a point cloud coloured by energy (drag to rotate, scroll to zoom).

```python exec="1" html="1"
import html
import math

import numpy as np
import pyvista as pv
from vtkmodules.vtkIOHDF import vtkHDFReader

import yamc

pv.OFF_SCREEN = True
source = yamc.NeutronSource(
    position=yamc.sources.CylindricalRing(
        radius=yamc.sources.Normal(300, 25),        # fuzzy major radius (cm)
        phi=yamc.sources.Uniform(0, 2 * math.pi),   # full ring
        z=yamc.sources.Normal(0, 20),               # spread about the midplane
    ),
    energy=yamc.sources.Normal(14.06e6, 3e5),
)
source.to_vtkhdf("source_points.vtkhdf", samples=4000)

reader = vtkHDFReader()
reader.SetFileName("source_points.vtkhdf")
reader.Update()
points = pv.wrap(reader.GetOutput())

plotter = pv.Plotter(window_size=(760, 560))
plotter.add_mesh(points, scalars="energy", cmap="plasma", render_points_as_spheres=True,
                 point_size=5, show_scalar_bar=False)
plotter.add_axes()
plotter.camera_position = "iso"

doc = plotter.export_html(filename=None).read()
print(
    f'<iframe srcdoc="{html.escape(doc, quote=True)}" width="100%" height="580" '
    'style="border:1px solid #ddd;border-radius:4px"></iframe>'
)
```



