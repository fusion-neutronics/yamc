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
# Fixed point -- the (x, y, z) tuple shorthand
source = yamc.NeutronSource(position=(0, 0, 0), energy=14.06e6)

# Cylindrical ring source (e.g. tokamak plasma)
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

### Energy distributions

```python
# Monoenergetic
energy = 14.06e6

# Discrete line spectrum -- e.g. Co-60's two main decay gammas,
# 1.17 MeV and 1.33 MeV, emitted in ~equal intensity
energy = yamc.sources.Discrete([1.1732e6, 1.3325e6], [0.5, 0.5])

# Histogram (piecewise-constant) spectrum: n+1 ascending bin edges (eV) and one
# per-bin weight. Weights are per-bin probability mass and are normalized
# internally, so a raw multigroup flux per group can be passed directly. Here
# 30% of neutrons are born uniformly in [0, 1] MeV and 70% in [1, 20] MeV.
energy = yamc.sources.Histogram([0.0, 1e6, 20e6], [0.3, 0.7])

# The bin edges can also be a built-in group structure name instead of a list
# (e.g. "CCFE-709", "VITAMIN-J-175"); pass one weight per group.
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
direction = yamc.sources.Isotropic()                 # default -- usually omitted
direction = yamc.sources.Monodirectional([0, 0, 1])  # beam in +z
```

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



