import constructive_solid_geometry_for_mc as csg4mc

def test_region_contains_debug():
    s = csg4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1)
    region = -s
    cell = csg4mc.Cell(cell_id=1, region=region)
    # Should be inside
    assert cell.contains(0.0, 0.0, 0.0)
    # Should be outside
    assert not cell.contains(3.0, 0.0, 0.0)
    # Should be on surface (may be inside or outside depending on convention)
    print('contains(2,0,0):', cell.contains(2.0, 0.0, 0.0))
    # Print region type for debugging
    print('region:', region)
