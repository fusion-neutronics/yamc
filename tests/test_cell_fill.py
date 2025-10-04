import materials_for_mc as csg

def test_cell_fill():
    surf1 = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=1)
    mat = csg.Material("fuel")
    mat.add_nuclide('Li6', 1.0)
    mat.set_density('g/cm3', 1.0)
    mat.read_nuclides_from_json("tendl-21")
    cell = csg.Cell(cell_id=1, region=-surf1, fill=mat)
    cell.fill.macroscopic_cross_section(reaction=1)
    assert cell.fill is not None
    assert cell.fill.name == "fuel"

def test_cell_fill_optional():
    surf1 = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=1)
    cell = csg.Cell(cell_id=2, region=-surf1)
    assert cell.fill is None
