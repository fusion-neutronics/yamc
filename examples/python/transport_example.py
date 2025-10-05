import materials_for_mc as mc

# mc.Config.set_cross_sections("tendl-21")
# Create sphere surface with vacuum boundary
sphere1 = mc.Sphere(
    surface_id=1,
    x0=0.0,
    y0=0.0,
    z0=0.0,
    r=1.0,
    boundary_type='Vacuum',
)
sphere2 = mc.Sphere(
    surface_id=1,
    x0=0.0,
    y0=0.0,
    z0=0.0,
    r=2.0,
    boundary_type='Vacuum',
)
region1 = -sphere1
region2 = -sphere2 & +sphere1

material1 = mc.Material()
material1.add_nuclide("Li6", 1.0)
material1.set_density("g/cm3", 5.5)
material1.read_nuclides_from_json({"Li6": "tests/Li6.json"})

material2 = mc.Material()
material2.add_nuclide("Be9", 1.0)
material2.set_density("g/cm3", 2.0)
material2.read_nuclides_from_json({"Be9": "tests/Be9.json"})

cell1 = mc.Cell(
    cell_id=1,
    name="sphere_cell",
    region=region1,
    fill=material1,
)
cell2 = mc.Cell(
    cell_id=1,
    name="sphere_cell",
    region=region1,
    fill=material2,
)
geometry = mc.Geometry(cells=[cell1, cell2])

source = mc.Source([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1e6)
settings = mc.Settings(particles=5, batches=3, source=source)

filter1 = mc.CellFilter(cell1)
tally1= mc.Tally()
tally1.filters = [filter1]
tally1.score=101
tally1.name="absorption in cell 1"

tally2= mc.Tally()
tally2.score=101
tally2.name="absorption in whole model"

tallies = [tally1, tally2]


model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
leakage_tally, tally1, tally2 = model.run()

# print(leakage_tally)
print(tally1, end="\n\n")
print(tally2)