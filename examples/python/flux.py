import materials_for_mc as mc
import time

# Create spherical geometry with Li6
sphere1 = mc.Surface(
    surface_id=1,
    kind=mc.SurfaceKind.sphere(0.0, 0.0, 0.0, 1.0),
    boundary_type=mc.BoundaryType.Transmission,
)

sphere2 = mc.Surface(
    surface_id=2,
    kind=mc.SurfaceKind.sphere(0.0, 0.0, 0.0, 200.0),
    boundary_type=mc.BoundaryType.Vacuum,
)

region1 = mc.Region.from_halfspace(mc.HalfspaceType.below(sphere1))
region2_inner = mc.Region.from_halfspace(mc.HalfspaceType.above(sphere1))
region2_outer = mc.Region.from_halfspace(mc.HalfspaceType.below(sphere2))
region2 = region2_inner.intersection(region2_outer)

# Create Li6 material
material = mc.Material()
material.set_density("g/cm3", 2.0)
material.add_nuclide("Li6", 1.0)
nuclide_json_map = {"Li6": "tests/Li6.json"}
material.read_nuclides_from_json(nuclide_json_map)

# Create cells
cell1 = mc.Cell(
    cell_id=1,
    region=region1,
    name="inner_sphere",
    material=None,
)

cell2 = mc.Cell(
    cell_id=2,
    region=region2,
    name="outer_annular",
    material=material,
)

geometry = mc.Geometry(cells=[cell1, cell2])

# Create source
source = mc.IndependentSource(
    space=[0.0, 0.0, 0.0],
    angle=mc.AngularDistribution.isotropic(),
    energy=14.06e6,
)

settings = mc.Settings(
    particles=10000,
    batches=10,
    source=source,
    seed=1,
)

# Create tallies: one for flux, one for tritium production
cell_filter = mc.CellFilter.new(cell2)

flux_tally = mc.Tally()
flux_tally.add_filter(cell_filter)
flux_tally.add_score(mc.FLUX_SCORE)  # Use FLUX_SCORE constant
flux_tally.set_name("flux")

tbr_tally = mc.Tally()
tbr_tally.add_filter(cell_filter)
tbr_tally.add_score(105)  # MT 105: n,t (tritium production)
tbr_tally.set_name("tbr")

# Create and run model
model = mc.Model(
    geometry=geometry,
    settings=settings,
    tallies=[flux_tally, tbr_tally],
)

start = time.time()
model.run()
elapsed = time.time() - start

print(f"\nSimulation completed in {elapsed:.6f} seconds.")

# Print results
flux_mean = flux_tally.get_mean()
tbr_mean = tbr_tally.get_mean()

print(f"Flux (track-length): {flux_mean[0]:.6e} cm")
print(f"TBR (tritium breeding ratio): {tbr_mean[0]:.6f}")
