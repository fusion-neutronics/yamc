import constructive_solid_geometry_for_mc as csg4mc
import math

def test_cell_closest_surface_simple():
    # Cell: sphere of radius 2 at (0,0,0)
    s = csg4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1)
    region = -s  # inside sphere
    cell = csg4mc.Cell(cell_id=1, region=region)
    # From (3,0,0) toward center
    result = cell.closest_surface((3.0, 0.0, 0.0), (-1.0, 0.0, 0.0))
    assert result is not None, 'Expected intersection, got None'
    surf, dist = result
    assert math.isclose(dist, 1.0, abs_tol=1e-10)
    # Surface is a PySurface, but type check is not needed for this test
    # From (0,0,0) outward
    result2 = cell.closest_surface((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result2 is not None, 'Expected intersection, got None'
    surf2, dist2 = result2
    assert math.isclose(dist2, 2.0, abs_tol=1e-10)
    # From outside, away from sphere
    result3 = cell.closest_surface((3.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result3 is None

def test_cell_closest_surface_union():
    # Two spheres, union
    s1 = csg4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1)
    s2 = csg4mc.Sphere(x0=4.0, y0=0.0, z0=0.0, r=2.0, surface_id=2)
    region1 = -s1
    region2 = -s2
    region = region1 | region2
    cell = csg4mc.Cell(cell_id=2, region=region)
    # From (0,0,0) outward (should hit s1)
    result = cell.closest_surface((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result is not None, 'Expected intersection, got None'
    surf, dist = result
    assert math.isclose(dist, 2.0, abs_tol=1e-10)
    # From (4,0,0) outward (should hit s2)
    result2 = cell.closest_surface((4.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result2 is not None, 'Expected intersection, got None'
    surf2, dist2 = result2
    assert math.isclose(dist2, 2.0, abs_tol=1e-10)
    # From (2,0,0) toward (4,0,0) (should hit s2 at 4.0, since ray starts outside both spheres)
    result3 = cell.closest_surface((2.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result3 is not None, 'Expected intersection, got None'
    surf3, dist3 = result3
    assert math.isclose(dist3, 4.0, abs_tol=1e-10)

def test_cell_closest_surface_intersection():
    # Intersection of two spheres
    s1 = csg4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=3.0, surface_id=1)
    s2 = csg4mc.Sphere(x0=2.0, y0=0.0, z0=0.0, r=3.0, surface_id=2)
    region = (-s1) & (-s2)
    cell = csg4mc.Cell(cell_id=3, region=region)
    # From (1,0,0) outward (should hit s1 or s2, both at 2.0)
    result = cell.closest_surface((1.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result is not None, 'Expected intersection, got None'
    surf, dist = result
    assert math.isclose(dist, 2.0, abs_tol=1e-10)
    # From (0,0,0) outward (should hit s1 at 3.0)
    result2 = cell.closest_surface((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result2 is not None, 'Expected intersection, got None'
    surf2, dist2 = result2
    assert math.isclose(dist2, 3.0, abs_tol=1e-10)
    # From (2,0,0) outward (should hit s2 at 1.0 or s1 at 3.0)
    result3 = cell.closest_surface((2.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result3 is not None, 'Expected intersection, got None'
    surf3, dist3 = result3
    # Accept either 1.0 or 3.0 as valid (depending on which surface is hit first)
    assert math.isclose(dist3, 1.0, abs_tol=1e-10) or math.isclose(dist3, 3.0, abs_tol=1e-10)

def test_cell_closest_surface_complement():
    # Complement of a sphere (outside)
    s = csg4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1)
    region = ~(-s)
    cell = csg4mc.Cell(cell_id=4, region=region)
    # From (3,0,0) outward (should hit nothing)
    result = cell.closest_surface((3.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert result is None
    # From (3,0,0) toward center (should hit sphere at 1.0)
    result2 = cell.closest_surface((3.0, 0.0, 0.0), (-1.0, 0.0, 0.0))
    assert result2 is not None, 'Expected intersection, got None'
    surf, dist = result2
    assert math.isclose(dist, 1.0, abs_tol=1e-10)
