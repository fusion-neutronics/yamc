import materials_for_mc as mc

# mc.Config.set_cross_sections("tendl-21")
# Create sphere surface with vacuum boundary
sphere = mc.Sphere(
    surface_id=1,
    x0=0.0,
    y0=0.0,
    z0=0.0,
    r=2.0,
    boundary_type='Vacuum',
)
region = -sphere
material = mc.Material()
material.add_nuclide("Li6", 1.0)
material.set_density("g/cm3", 5.5)
# material.read_nuclides_from_json("tendl-21")
material.read_nuclides_from_json({"Li6": "tests/Li6.json"})
cell = mc.Cell(
    cell_id=1,
    name="sphere_cell",
    region=region,
    fill=material,
)
geometry = mc.Geometry(cells=[cell])

source = mc.Source([0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 1e6)
settings = mc.Settings(particles=2, batches=1, source=source)


tally1= mc.Tally()
tally1.score=101
tally1.name="absorption in whole model"
tallies = [tally1]
    

model = mc.Model(geometry=geometry, settings=settings, tally=tallies)
leakage_tally, absorption_tally = model.run()

print(leakage_tally)
print(absorption_tally)