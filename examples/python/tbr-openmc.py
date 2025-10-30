import materials_for_mc as mc
import openmc as mc
mc.config['cross_sections'] = "/home/jon/nuclear_data/cross_sections.xml"
# material1.read_nuclides_from_json({"Li6": "tests/Li6.json", "Li7": "tests/Li7.json"})
    
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
material1.add_nuclide("Li6", 0.07)  # Li4SiO4
material1.add_nuclide("Li7", 0.93)
# material1.add_element("O", 4.0)
# material1.add_element("Si", 1.0)
material1.set_density("g/cm3", 2.0)
# material1.read_nuclides_from_json({"Li6": "tests/Li6.json", "Li7": "tests/Li7.json"})


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

source = mc.IndependentSource(space=mc.stats.Point([0,0,0]), angle=mc.stats.Isotropic(), energy=mc.stats.Discrete([14060000.0], [1.0]))
settings = mc.Settings(particles=500, batches=200, source=source, seed=1, run_mode='fixed source')  # Increased for better statistics

# Create tallies with CellFilters
cell_filter2 = mc.CellFilter(cell2)
tally1 = mc.Tally()
tally1.filters = [cell_filter2]
tally1.scores = [205]  # n,t
tally1.name = "tbr"
tallies = [tally1]

model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
model.run(apply_tally_results=True)

# Tallies are updated in place!
print(f"TBR (tritium breeding ratio): {tally1.mean}")
