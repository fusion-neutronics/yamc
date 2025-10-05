import materials_for_mc as mc

# mc.Config.set_cross_sections("tendl-21")
# Create sphere surface with transmission boundary 
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
    r=2.0,
    boundary_type='vacuum',
)
region1 = -sphere1
region2 = +sphere1 & -sphere2

material1 = mc.Material()
material1.add_nuclide("Li6", 1.0)
material1.set_density("g/cm3", 10.0)  # Higher density for more absorption
material1.read_nuclides_from_json({"Li6": "tests/Li6.json"})

material2 = mc.Material()
material2.add_nuclide("Be9", 1.0)
material2.set_density("g/cm3", 20.0)  # Higher density for more absorption
material2.read_nuclides_from_json({"Be9": "tests/Be9.json"})

cell1 = mc.Cell(
    cell_id=1,
    name="sphere_cell",
    region=region1,
    fill=material1,
)
cell2 = mc.Cell(
    cell_id=2,
    name="annular_cell",
    region=region2,
    fill=material2,
)
geometry = mc.Geometry(cells=[cell1, cell2])

source = mc.Source(position=[0.0, 0.0, 0.0], direction=[0.0, 0.0, 1.0], energy=1e6)
settings = mc.Settings(particles=100, batches=10, source=source)

filter1 = mc.CellFilter(cell1)
tally1= mc.Tally()
tally1.filters = [filter1]
tally1.score=101
tally1.name="absorption in cell 1"

filter2 = mc.CellFilter(cell2)
tally2= mc.Tally()
tally2.filters = [filter2]
tally2.score=101
tally2.name="absorption in cell 2"

tally3= mc.Tally()
tally3.score=101
tally3.name="absorption in whole model"

tallies = [tally1, tally2, tally3]


model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
leakage_tally, tally1, tally2, tally3 = model.run()

print(leakage_tally)
print(tally1, end="\n\n")
print(tally2, end="\n\n")
print(tally3, end="\n\n")

assert tally1.mean != tally2.mean
assert abs((tally1.mean + tally2.mean) - tally3.mean) < 1e-10 # checking sum of absorptions in cells equals total absorption
assert abs((1 - tally3.mean) - leakage_tally.mean) < 1e-10 # checking particles are either absorbed or leaked