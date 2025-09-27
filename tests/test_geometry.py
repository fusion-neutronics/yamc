import constructive_solid_geometry_for_mc as csg

def test_geometry_find_cell():
    surf = csg.Sphere(x0=0, y0=0, z0=0, r=2, surface_id=1)
    region = -surf
    cell = csg.Cell(cell_id=1, region=region)
    geometry = csg.Geometry([cell])
    assert geometry.find_cell(0, 0, 0) is not None
    assert geometry.find_cell(5, 0, 0) is None
