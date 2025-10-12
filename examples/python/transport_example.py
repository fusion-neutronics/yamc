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
material2.add_nuclide("Li7", 1.0)
material2.set_density("g/cm3", 20.0)  # Higher density for more absorption
material2.read_nuclides_from_json({"Li7": "tests/Li7.json"})

cell1 = mc.Cell(
    cell_id=1,
    # name="sphere_cell",
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

source = mc.IndependentSource()  # defaults to 14.06MeV energy, space = 0,0,0 direction = isotropic
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

tally4= mc.Tally()
tally4.score=4
tally4.name="inelastic in whole model"

tally5= mc.Tally()
tally5.score=55
tally5.name="inelastic (55) in whole model"

tallies = [tally1, tally2, tally3, tally4, tally5]

# Add fission tally (MT=18)
fission_tally = mc.Tally()
fission_tally.score = 18
fission_tally.name = "fission in whole model"
tallies.append(fission_tally)


model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)

# Unpack tallies including fission tally
leakage_tally, tally1, tally2, tally3, tally4, tally5, fission_tally = model.run()

print(leakage_tally)
print(tally1, end="\n\n")
print(tally2, end="\n\n")
print(tally3, end="\n\n")
print(tally4, end="\n\n")
print(tally5, end="\n\n")
print(fission_tally, end="\n\n")

assert tally1.mean != tally2.mean
assert abs((tally1.mean + tally2.mean) - tally3.mean) < 1e-10 # checking sum of absorptions in cells equals total absorption
assert fission_tally.mean == 0.0  # No fission should occur in Li6 or Li7
assert abs((1 - tally3.mean - fission_tally.mean) - leakage_tally.mean) < 1e-10 # checking particles are either absorbed, fissioned, or leaked