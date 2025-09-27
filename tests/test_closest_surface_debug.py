import constructive_solid_geometry_for_mc as csg4mc
import math

def test_closest_surface_debug():
    s = csg4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1)
    region = -s
    cell = csg4mc.Cell(cell_id=1, region=region)
    # From (3,0,0) toward center
    result = cell.closest_surface((3.0, 0.0, 0.0), (-1.0, 0.0, 0.0))
    print('closest_surface((3,0,0),(-1,0,0)):', result)
    # From (0,0,0) outward
    result2 = cell.closest_surface((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    print('closest_surface((0,0,0),(1,0,0)):', result2)
    # From (3,0,0) outward (should be None)
    result3 = cell.closest_surface((3.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    print('closest_surface((3,0,0),(1,0,0)):', result3)
    # From (2,0,0) outward (on surface)
    result4 = cell.closest_surface((2.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    print('closest_surface((2,0,0),(1,0,0)):', result4)
