import pytest
from constructive_solid_geometry_for_mc import XPlane, YPlane, ZPlane, Sphere, Cylinder, ZCylinder

def test_xplane_creation():
    s = XPlane(x0=1.0, surface_id=42)
    assert s.id == 42
    assert s.evaluate((1.0, 0.0, 0.0)) == pytest.approx(0.0)

def test_yplane_creation():
    s = YPlane(y0=2.0, surface_id=43)
    assert s.id == 43
    assert s.evaluate((0.0, 2.0, 0.0)) == pytest.approx(0.0)

def test_zplane_creation():
    s = ZPlane(z0=3.0, surface_id=44)
    assert s.id == 44
    assert s.evaluate((0.0, 0.0, 3.0)) == pytest.approx(0.0)

def test_sphere_creation():
    s = Sphere(x0=1.0, y0=2.0, z0=3.0, r=5.0, surface_id=45)
    assert s.id == 45
    assert s.evaluate((1.0, 2.0, 8.0)) == pytest.approx(0.0)

def test_cylinder_creation():
    s = Cylinder(x0=1.0, y0=2.0, z0=3.0, axis_x=0.0, axis_y=1.0, axis_z=0.0, r=2.0, surface_id=46)
    assert s.id == 46
    # Point at radius from origin, perpendicular to axis (Y axis)
    assert s.evaluate((3.0, 2.0, 3.0)) == pytest.approx(0.0)

def test_zcylinder_creation():
    s = ZCylinder(x0=1.0, y0=2.0, r=3.0, surface_id=47)
    assert s.id == 47
    # Point at radius from center in XY plane
    assert s.evaluate((4.0, 2.0, 0.0)) == pytest.approx(0.0)

def test_boundary_type_default():
    s = XPlane(x0=1.0, surface_id=42)
    assert str(s.boundary_type) == "vacuum"

def test_boundary_type_vacuum():
    s = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=1, boundary_type="transmission")
    assert str(s.boundary_type) == "transmission"

def test_set_boundary_type():
    s = ZCylinder(x0=0.0, y0=0.0, r=1.0, surface_id=2)
    assert str(s.boundary_type) == "vacuum"
    
    s.boundary_type = "transmission"
    assert str(s.boundary_type) == "transmission"
    
    s.boundary_type = "vacuum"
    assert str(s.boundary_type) == "vacuum"

def test_invalid_boundary_type():
    with pytest.raises(ValueError):
        s = XPlane(x0=1.0, surface_id=42, boundary_type="invalid")
        
def test_invalid_set_boundary_type():
    s = XPlane(x0=1.0, surface_id=42)
    with pytest.raises(ValueError):
        s.boundary_type = "invalid"
