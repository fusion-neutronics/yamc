import constructive_solid_geometry_for_mc as csg

def test_cell_fill():
    surf1 = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=1)
    mat = csg.Material("fuel")
    cell = csg.Cell(cell_id=1, region=-surf1, fill=mat)
    assert cell.fill is not None
    assert cell.fill.name == "fuel"

def test_cell_fill_optional():
    surf1 = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=1)
    cell = csg.Cell(cell_id=2, region=-surf1)
    assert cell.fill is None
