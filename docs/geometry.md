# Geometry

## CSG Geometry

YAMC uses constructive solid geometry (CSG) to define simulation domains. Surfaces partition space into half-spaces, which are combined with boolean operations to form regions, cells, and ultimately a geometry.

### Surfaces

All surfaces take an optional ``surface_id`` and ``boundary`` (``"transmission"`` or ``"vacuum"``).

**Planes** are infinite planar surfaces:

```python
# Axis-aligned planes
xp = yamc.Plane(axis="x", offset=5.0)
yp = yamc.Plane(axis="y", offset=0.0)
zp = yamc.Plane(axis="z", offset=-10.0)

# General plane: ax + by + cz = d
p = yamc.Plane(axis=(1.0, 1.0, 0.0), offset=5.0 / ((1.0**2 + 1.0**2 + 0.0**2) ** 0.5))
```

A **sphere** is centered at ``(x0, y0, z0)`` with radius ``r``:

```python
sphere = yamc.Sphere(x0=0.0, y0=0.0, z0=0.0, radius=10.0)
```

**Cylinders** are infinite along their axis:

```python
# Axis-aligned cylinders
zc = yamc.Cylinder(axis="z", x0=0.0, y0=0.0, radius=5.0)
xc = yamc.Cylinder(axis="x", y0=0.0, z0=0.0, radius=3.0)
yc = yamc.Cylinder(axis="y", x0=0.0, z0=0.0, radius=3.0)

# General cylinder (arbitrary axis)
cyl = yamc.Cylinder(axis=(1.0, 1.0, 0.0), x0=0.0, y0=0.0, z0=0.0, radius=5.0)
```

For **tori**, ``r_major`` is the major radius and ``r_minor`` the minor
radius along the symmetry axis; pass ``r_minor_2`` for a different
transverse minor radius (an elliptical cross-section). ``axis`` is
``'x'``, ``'y'``, or ``'z'`` (arbitrary directions are not yet
supported):

```python
zt = yamc.Torus(axis="z", x0=0.0, y0=0.0, z0=0.0, r_major=10.0, r_minor=3.0)
xt = yamc.Torus(axis="x", r_major=10.0, r_minor=3.0)
yt = yamc.Torus(axis="y", r_major=10.0, r_minor=3.0)
```

**Cones** are double-sheeted, with the apex at ``(x0, y0, z0)`` and
``opening_angle`` the half-opening angle in degrees. Cap with planes
to keep one sheet:

```python
import math

# Axis-aligned cones
zc = yamc.Cone(axis="z", x0=0.0, y0=0.0, z0=0.0, opening_angle=math.degrees(math.atan(0.25**0.5)))
xc = yamc.Cone(axis="x", opening_angle=math.degrees(math.atan(1.0**0.5)))
yc = yamc.Cone(axis="y", opening_angle=math.degrees(math.atan(1.0**0.5)))

# General cone (arbitrary axis)
cone = yamc.Cone(axis=(0.0, 1.0, 0.0), x0=0.0, y0=0.0, z0=0.0, opening_angle=math.degrees(math.atan(0.5**0.5)))
```

The **general quadric** is a catch-all primitive,
``a x^2 + b y^2 + c z^2 + d xy + e yz + f xz + g x + h y + j z + k = 0``.
The coefficient names skip ``i``. It is useful for porting existing
models and for surfaces with no dedicated variant. The dedicated
variants are faster and give finite bounding boxes, so prefer
them when one fits:

```python
# A sphere of radius 10 written as a quadric
q = yamc.Quadric(a=1.0, b=1.0, c=1.0, k=-100.0)
```

**Boundary types** are set on the outermost surfaces:

- ``"transmission"`` (default): particles pass through freely
- ``"vacuum"``: particles are killed on crossing

```python
outer = yamc.Sphere(radius=50.0, boundary="vacuum")
```

### Halfspaces and regions

Surfaces divide space into two halfspaces. The ``.below`` property gives the inside/below side (negative sense), ``.above`` gives the outside/above side (positive sense).

```python
sphere = yamc.Sphere(radius=10.0)
inside = sphere.below   # r < 10
outside = sphere.above  # r > 10
```

### Boolean operators

Regions are combined with ``&`` (intersection), ``|`` (union), and ``~`` (complement) to build more complicated shapes.

**Intersection** (``&``) is the space where both regions overlap:

```python
# A cylinder capped by two planes gives a finite cylinder
cylinder = yamc.Cylinder(axis="z", radius=5.0)
bottom = yamc.Plane(axis="z", offset=-10.0)
top = yamc.Plane(axis="z", offset=10.0)
capped_cylinder = cylinder.below & bottom.above & top.below
```

**Union** (``|``) is the space belonging to either region:

```python
# Two cylinders forming a cross
cyl_x = yamc.Cylinder(axis="x", radius=2.0)
cyl_z = yamc.Cylinder(axis="z", radius=2.0)
cross = cyl_x.below | cyl_z.below
```

**Complement** (``~``) is everything outside a region:

```python
# Hollow cylinder (pipe)
inner_cyl = yamc.Cylinder(axis="z", radius=3.0)
outer_cyl = yamc.Cylinder(axis="z", radius=5.0)
bottom = yamc.Plane(axis="z", offset=-10.0)
top = yamc.Plane(axis="z", offset=10.0)
pipe = ~inner_cyl.below & outer_cyl.below & bottom.above & top.below
```

**Combining them** gives, for example, a box with a cylindrical hole:

```python
# Box
xmin, xmax = yamc.Plane(axis="x", offset=-8.0), yamc.Plane(axis="x", offset=8.0)
ymin, ymax = yamc.Plane(axis="y", offset=-8.0), yamc.Plane(axis="y", offset=8.0)
zmin, zmax = yamc.Plane(axis="z", offset=-8.0), yamc.Plane(axis="z", offset=8.0)
box = xmin.above & xmax.below & ymin.above & ymax.below & zmin.above & zmax.below

# Cylindrical hole through the box
hole = yamc.Cylinder(axis="z", radius=3.0)
box_with_hole = box & ~hole.below  # box minus cylinder interior
```

### Cells

A cell combines a region with an optional material fill. Cells with no fill are void.

<!-- doctest: skip -->
```python
steel = yamc.Material(composition={"Fe": 1.0}, density=7.874, name="steel")
water = yamc.Material(composition={"H2O": 1.0}, density=1.0, name="water")

# Steel pipe filled with water, inside a vacuum box
pipe_cell = yamc.Cell(region=pipe, material=steel, name="pipe")
water_cell = yamc.Cell(region=inner_cyl.below & bottom.above & top.below, material=water, name="coolant")
void_cell = yamc.Cell(region=box & outer_cyl.above, name="void")
```

### Geometry

A ``Geometry`` is constructed from a list of cells. Cell IDs are auto-assigned if not set. Material IDs are not: every ``Material`` has one from when it was built, and a ``Geometry`` refuses two different materials with the same ID.

<!-- doctest: skip -->
```python
geometry = yamc.Geometry([pipe_cell, water_cell, void_cell])
```

### Plotting

The interactive plot renders in Jupyter notebooks and can be saved as HTML (with pan/zoom/slice controls) or PNG. ``Cell`` and ``Region`` objects also support ``.plot()``.

<!-- doctest: skip -->
```python
plot = geometry.plot(basis="xz", color_by="material", outline="cell", resolution=400000)
plot.save("geometry.html")   # interactive HTML
plot.save("geometry.png")    # static image

# Also works on individual cells and regions
plot = cell.plot(basis="xz", resolution=400000)
plot = region.plot(basis="xz", resolution=400000)
```

The geometry viewer below is rendered live when the docs are built (drag to pan, scroll to zoom):

```python exec="1" html="1"
import html
import yamc

steel = yamc.Material(composition={"Fe": 1.0}, density=7.874, name="steel")
water = yamc.Material(composition={"H2O": 1.0}, density=1.0, name="water")
inner, outer = yamc.Cylinder(axis="z", radius=2.0), yamc.Cylinder(axis="z", radius=8.0)
bottom, top = yamc.Plane(axis="z", offset=-12.0), yamc.Plane(axis="z", offset=12.0)
box = (
    yamc.Plane(axis="x", offset=-15.0, boundary="vacuum").above
    & yamc.Plane(axis="x", offset=15.0, boundary="vacuum").below
    & yamc.Plane(axis="y", offset=-15.0, boundary="vacuum").above
    & yamc.Plane(axis="y", offset=15.0, boundary="vacuum").below
    & yamc.Plane(axis="z", offset=-15.0, boundary="vacuum").above
    & yamc.Plane(axis="z", offset=15.0, boundary="vacuum").below
)
pipe = yamc.Cell(region=inner.above & outer.below & bottom.above & top.below, material=steel, name="pipe")
coolant = yamc.Cell(region=inner.below & bottom.above & top.below, material=water, name="coolant")
void = yamc.Cell(region=box & ~(outer.below & bottom.above & top.below), name="void")
geometry = yamc.Geometry([pipe, coolant, void])

doc = geometry.plot(basis="xz", color_by="material", outline="cell")._repr_html_()
print(f'<iframe srcdoc="{html.escape(doc, quote=True)}" width="100%" height="520" style="border:1px solid #ddd;border-radius:4px"></iframe>')
```


For manual control with matplotlib, use ``sample_slice()`` (the figure below is rendered live when the docs are built):

```python exec="1" html="1" source="above"
import numpy as np
import matplotlib.pyplot as plt
import yamc

# Steel pipe + water coolant inside a vacuum box
steel = yamc.Material(composition={"Fe": 1.0}, density=7.874, temperature=294, name="steel")
water = yamc.Material(composition={"H2O": 1.0}, density=1.0, temperature=294, name="water")
inner_cyl, outer_cyl = yamc.Cylinder(axis="z", radius=2.0), yamc.Cylinder(axis="z", radius=8.0)
bottom, top = yamc.Plane(axis="z", offset=-12.0), yamc.Plane(axis="z", offset=12.0)
box = (yamc.Plane(axis="x", offset=-15.0, boundary="vacuum").above & yamc.Plane(axis="x", offset=15.0, boundary="vacuum").below
       & yamc.Plane(axis="y", offset=-15.0, boundary="vacuum").above & yamc.Plane(axis="y", offset=15.0, boundary="vacuum").below
       & yamc.Plane(axis="z", offset=-15.0, boundary="vacuum").above & yamc.Plane(axis="z", offset=15.0, boundary="vacuum").below)
pipe = yamc.Cell(region=inner_cyl.above & outer_cyl.below & bottom.above & top.below, material=steel, name="pipe")
coolant = yamc.Cell(region=inner_cyl.below & bottom.above & top.below, material=water, name="coolant")
void = yamc.Cell(region=box & ~(outer_cyl.below & bottom.above & top.below), name="void")
geometry = yamc.Geometry([pipe, coolant, void])

slice_data = geometry.sample_slice(basis="xz", resolution=(400, 400))
cell_ids = np.array(slice_data.cell_ids).reshape(400, 400)

fig, ax = plt.subplots()
ax.imshow(cell_ids, extent=slice_data.extent, origin="lower", cmap="Set2")
ax.set_xlabel("x [cm]")
ax.set_ylabel("z [cm]")
ax.set_title("CSG geometry slice (xz)")

from io import StringIO as _S
_b = _S(); plt.savefig(_b, format="svg", bbox_inches="tight"); print(_b.getvalue())
```

### VTK-HDF geometry export

Geometry can be exported to VTK-HDF ImageData files for visualization in [ParaView](https://www.paraview.org/). ``Cell`` and ``Region`` objects also support ``to_vtkhdf()``.

<!-- doctest: skip -->
```python
geometry.to_vtkhdf("geometry.vtkhdf", resolution=10000)
cell.to_vtkhdf("cell.vtkhdf")
region.to_vtkhdf("region.vtkhdf")
```

The viewer below is rendered live when these docs are built: it voxelizes the steel-pipe-and-water geometry to ``geometry.vtkhdf``, reads it back, and shows it cut in half so the material regions are visible, coloured by ``material_id`` (drag to rotate, scroll to zoom).

```python exec="1" html="1"
import html

import pyvista as pv
from vtkmodules.vtkIOHDF import vtkHDFReader

import yamc

pv.OFF_SCREEN = True
steel = yamc.Material(composition={"Fe": 1.0}, density=7.874, name="steel")
water = yamc.Material(composition={"H2O": 1.0}, density=1.0, name="water")
inner, outer = yamc.Cylinder(axis="z", radius=2.0), yamc.Cylinder(axis="z", radius=8.0)
bottom, top = yamc.Plane(axis="z", offset=-12.0), yamc.Plane(axis="z", offset=12.0)
box = (
    yamc.Plane(axis="x", offset=-15.0, boundary="vacuum").above
    & yamc.Plane(axis="x", offset=15.0, boundary="vacuum").below
    & yamc.Plane(axis="y", offset=-15.0, boundary="vacuum").above
    & yamc.Plane(axis="y", offset=15.0, boundary="vacuum").below
    & yamc.Plane(axis="z", offset=-15.0, boundary="vacuum").above
    & yamc.Plane(axis="z", offset=15.0, boundary="vacuum").below
)
pipe = yamc.Cell(region=inner.above & outer.below & bottom.above & top.below, material=steel, name="pipe")
coolant = yamc.Cell(region=inner.below & bottom.above & top.below, material=water, name="coolant")
void = yamc.Cell(region=box & ~(outer.below & bottom.above & top.below), name="void")
geometry = yamc.Geometry([pipe, coolant, void])
geometry.to_vtkhdf("geometry.vtkhdf", resolution=300000)

reader = vtkHDFReader()
reader.SetFileName("geometry.vtkhdf")
reader.Update()
grid = pv.wrap(reader.GetOutput())
solid = grid.threshold(-0.5, scalars="material_id").clip(normal="x")  # drop void, cut in half

plotter = pv.Plotter(window_size=(760, 560))
plotter.add_mesh(solid, scalars="material_id", cmap="Set2", show_scalar_bar=False)
plotter.add_axes()
plotter.camera_position = "iso"

doc = plotter.export_html(filename=None).read()
print(
    f'<iframe srcdoc="{html.escape(doc, quote=True)}" width="100%" height="580" '
    'style="border:1px solid #ddd;border-radius:4px"></iframe>'
)
```


### Volume calculation

Stochastic volume calculation samples random points and checks which cell
contains each one. Use it to find the volumes of CSG cells.

<!-- doctest: skip -->
```python
# Single cell
result = cell.calculate_volume(samples=1_000_000)
print(f"Volume: {result.volume:.2f} +/- {result.std_dev:.2f} cm^3")

# All cells in geometry
volumes = geometry.calculate_volume(samples=1_000_000)
for cell_id, result in volumes.items():
    print(f"Cell {cell_id}: {result.volume:.2f} cm^3")
```

Meshed geometries do not need stochastic volume calculation. Each volume's
size is computed analytically (via the divergence theorem over its bounding
triangles) when the mesh is loaded, and is available as
``MeshGeometry.volume_measures``.

<!-- doctest: skip -->
```python
mesh_geom = yamc.MeshGeometry(filename="geometry.arrow", materials={"mat:fuel": fuel})

for name, vol in zip(mesh_geom.volume_names, mesh_geom.volume_measures):
    print(f"{name}: {vol:.2f} cm^3")
```

For consistency with the CSG API, ``MeshGeometry.calculate_volume()`` is
also available. It returns the same dict-of-VolumeResult shape as
``Geometry.calculate_volume()`` but with exact values (``std_dev=0``)
computed from the analytic measures, so the same calling code works
for both CSG and mesh geometries.

<!-- doctest: skip -->
```python
volumes = mesh_geom.calculate_volume()  # samples= and bounding_box= are ignored
for vol_id, result in volumes.items():
    print(f"Volume {vol_id}: {result.volume:.2f} cm^3 (exact)")
```


## CAD Geometry

YAMC includes a CAD-to-mesh pipeline in the ``yamc.cad`` module (install with ``pip install "yamc[cad]"``). [CadQuery](https://cadquery.readthedocs.io/) is a Python library for building parametric 3D CAD models with a script-based API. The pipeline lets CadQuery geometry be used directly for neutron transport.

CadQuery can also load externally created STEP files (e.g. exported from another CAD tool) with ``cq.importers.importStep("part.step")``. The imported CAD should be clean and watertight so that it meshes correctly and works in particle transport. Avoid gaps, overlaps and self-intersections between solids; they cause meshing failures or lost or double-counted particles.

### Defining geometry in CadQuery

```python
import cadquery as cq

sphere = cq.Workplane("XY").sphere(10.0)
box = cq.Workplane("XY").box(30, 30, 30).cut(sphere)  # no overlap: every point in one solid

assembly = cq.Assembly()
assembly.add(sphere, name="fuel", color=cq.Color("red"))
assembly.add(box, name="moderator", color=cq.Color("blue"))
```

### Meshing (CadToYamc)

``CadToYamc`` imprints the assembly by default, meshes surfaces with the yamm scene mesher, and optionally fills volumes with tetrahedra using the yamm tet mesher.

<!-- doctest: skip -->
```python
from yamc.cad import CadToYamc

mesher = CadToYamc()
mesher.add_cadquery_object(assembly, material_tags=["fuel", "moderator"])

# Surface-only mesh
mesh = mesher.mesh(tolerance=0.1, angular_tolerance=0.1)

# Per-solid tolerances, keyed by material tag
mesh = mesher.mesh(tolerance={"fuel": 0.01, "moderator": 0.5})

# Selective tet meshing (specific volumes)
mesh = mesher.mesh(tet_volumes=["fuel"], target_edge_length=2.0)
```

The ``material_tags`` can also be ``"assembly_names"`` or ``"assembly_materials"`` to derive tags from the assembly structure.

``tolerance`` (chordal deflection) and ``angular_tolerance`` (radians) each take one number for every solid, or a dict keyed by material tag (every solid with that tag) or by 1-based solid ID (the keys of ``mesh.solid_faces``). Solids the dict does not name take 0.1. A surface two solids share is meshed at the finer of their two values, so the mesh stays watertight, and the ``mesh()`` docstring lists the other cases where a solid is meshed finer than it asks.

### Meshing without imprinting

Imprinting merges the surfaces where solids touch, so the mesh is conformal
across them. It is also the slow part of meshing a large assembly, and it
rejects CAD that will not imprint cleanly. ``imprint=False`` skips it and
meshes every solid on its own:

<!-- doctest: skip -->
```python
mesh = mesher.mesh(imprint=False, contact_tolerance=0.05)
```

Touching solids then meet along two separately faceted surfaces that may
overlap or leave a thin gap. ``MeshGeometry`` tracks such a mesh
discontinuously: a particle leaving one solid is placed in whichever touching
solid contains it, instead of being read off a shared surface. Solids closer
than ``contact_tolerance`` count as touching, so an overlap or a gap thinner
than it is crossed straight into the next solid, with an overlap staying with
the solid the particle was already in. A gap that is in the CAD stays a gap
whatever the tolerance: the pair is checked against the CAD and recorded as
separated. ``contact_tolerance`` defaults to twice ``tolerance`` (with
per-solid values, twice the largest), which covers the slivers that two
facetings of one CAD surface leave between them; raise it to accept sloppier
CAD.

Two solids overlapping deeper than ``contact_tolerance``, or without it twice
the larger ``tolerance`` of the two, is an error by
default. ``overlaps="smaller_wins"`` gives each such overlap to the smaller
solid of the pair, wherever a particle arrives from, and takes the shared
volume off the larger one. ``overlaps="ignore"`` skips the check.

The mesh file records how it was made, and ``MeshGeometry`` reports it:

<!-- doctest: skip -->
```python
print(mesh_geom.discontinuous)      # True for a mesh made without imprinting
print(mesh_geom.contact_tolerance)  # cm; None for an imprinted mesh
```

### Export

<!-- doctest: skip -->
```python
# Arrow IPC (for YAMC transport)
mesher.to_arrow("geometry.arrow")

# VTKHDF (for ParaView visualization)
mesher.to_vtkhdf("geometry.vtkhdf")
```

### MeshGeometry

``MeshGeometry`` loads a mesh file and assigns materials to each volume by name.

<!-- doctest: skip -->
```python
fuel = yamc.Material(composition={"U235": 0.05, "U238": 0.95}, density=19.1)
moderator = yamc.Material(composition={"H2O": 1.0}, density=1.0)

mesh_geom = yamc.MeshGeometry(
    filename="geometry.arrow",
    materials={
        "mat:fuel": fuel,
        "mat:moderator": moderator,
    },
)
```

The mesh file is a YAMC Arrow IPC mesh (``.arrow``); other meshers are supported by writing that format. Material keys use the ``"mat:<name>"`` convention matching the export tags.

Optional parameters control the vacuum boundary:

<!-- doctest: skip -->
```python
mesh_geom = yamc.MeshGeometry(
    filename="geometry.arrow",
    materials={"mat:fuel": fuel},
    implicit_complement_material=None,  # void outside mesh
    graveyard_offset=10.0,              # vacuum shell 10 cm outside mesh
)
```

``contact_tolerance`` overrides the distance the file records and tracks the mesh discontinuously, even if it was written imprinted.

A whole-model ``MeshGeometry`` with no vacuum boundary (no ``graveyard_offset``
and no surface with a vacuum condition) is closed before transport: with a void
implicit complement it gets a vacuum box 10 cm outside its bounding box, and with
a material in the implicit complement the run is refused, since where the box
sits would decide how much of that material surrounds the model. Pass
``graveyard_offset`` to choose the margin. A mesh used as the fill of a CSG cell
is bounded by the CSG geometry and gets nothing added.

Properties:

<!-- doctest: skip -->
```python
print(mesh_geom.num_volumes)      # number of mesh volumes
print(mesh_geom.volume_names)     # material name per volume
print(mesh_geom.volume_measures)  # analytical volume per volume [cm^3]
```

Plotting works the same as CSG geometry:

<!-- doctest: skip -->
```python
plot = mesh_geom.plot(basis="xz", color_by="material", outline="cell")
plot.save("mesh_geometry.html")
```

``MeshGeometry`` can be passed directly to ``Model()`` in place of a CSG ``Geometry``.


