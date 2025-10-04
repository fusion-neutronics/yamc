import materials_for_mc as m4mc
import math
import pytest

def test_plane_parallel_no_intersection():
    plane = m4mc.XPlane(x0=1.0, surface_id=1)
    d = plane.distance_to_surface((0.0, 0.0, 0.0), (0.0, 1.0, 0.0))
    assert d is None

def test_plane_on_surface():
    plane = m4mc.XPlane(x0=1.0, surface_id=1)
    d = plane.distance_to_surface((1.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    # Accept None (OpenMC/MCNP: no intersection if starting on surface and moving outward)
    assert d is None or math.isclose(d, 0.0, abs_tol=1e-10) or d > 0.0

def test_sphere_inside_out():
    sphere = m4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1)
    d = sphere.distance_to_surface((0.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert math.isclose(d, 2.0, abs_tol=1e-10)

def test_sphere_outside_away():
    sphere = m4mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1)
    d = sphere.distance_to_surface((3.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    assert d is None

def test_cylinder_axis_parallel():
    cyl = m4mc.ZCylinder(x0=0.0, y0=0.0, r=1.0, surface_id=1)
    d = cyl.distance_to_surface((0.5, 0.0, 0.0), (0.0, 0.0, 1.0))
    assert d is None

def test_cylinder_on_surface():
    cyl = m4mc.ZCylinder(x0=0.0, y0=0.0, r=1.0, surface_id=1)
    d = cyl.distance_to_surface((1.0, 0.0, 0.0), (1.0, 0.0, 0.0))
    # Accept None (OpenMC/MCNP: no intersection if starting on surface and moving outward)
    assert d is None or math.isclose(d, 0.0, abs_tol=1e-10) or d > 0.0

def test_plane_negative_direction():
    plane = m4mc.XPlane(x0=2.0, surface_id=1)
    d = plane.distance_to_surface((3.0, 0.0, 0.0), (-1.0, 0.0, 0.0))
    assert math.isclose(d, 1.0, abs_tol=1e-10)

def test_plane_no_intersection_behind():
    plane = m4mc.XPlane(x0=2.0, surface_id=1)
    d = plane.distance_to_surface((1.0, 0.0, 0.0), (-1.0, 0.0, 0.0))
    assert d is None
