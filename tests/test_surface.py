import pytest
from materials_for_mc import Plane, Geometry, Cell, Material, XPlane, YPlane, ZPlane, Sphere, Cylinder, ZCylinder

def test_xplane_creation():
    s = XPlane(x0=1.0, surface_id=42)
    assert s.surface_id == 42
    assert s.evaluate((1.0, 0.0, 0.0)) == pytest.approx(0.0)

def test_yplane_creation():
    s = YPlane(y0=2.0, surface_id=43)
    assert s.surface_id == 43
    assert s.evaluate((0.0, 2.0, 0.0)) == pytest.approx(0.0)

def test_zplane_creation():
    s = ZPlane(z0=3.0, surface_id=44)
    assert s.surface_id == 44
    assert s.evaluate((0.0, 0.0, 3.0)) == pytest.approx(0.0)

def test_sphere_creation():
    s = Sphere(x0=1.0, y0=2.0, z0=3.0, r=5.0, surface_id=45)
    assert s.surface_id == 45
    assert s.evaluate((1.0, 2.0, 8.0)) == pytest.approx(0.0)

def test_cylinder_creation():
    s = Cylinder(x0=1.0, y0=2.0, z0=3.0, axis_x=0.0, axis_y=1.0, axis_z=0.0, r=2.0, surface_id=46)
    assert s.surface_id == 46
    # Point at radius from origin, perpendicular to axis (Y axis)
    assert s.evaluate((3.0, 2.0, 3.0)) == pytest.approx(0.0)

def test_zcylinder_creation():
    s = ZCylinder(x0=1.0, y0=2.0, r=3.0, surface_id=47)
    assert s.surface_id == 47
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


def test_surface_optional_id():
    """Test that surface IDs are optional and default to None."""
    # Create surfaces without specifying surface_id
    sphere = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0)
    plane = Plane(a=1.0, b=0.0, c=0.0, d=0.0)
    cylinder = ZCylinder(x0=0.0, y0=0.0, r=1.0)
    
    # Surface IDs should default to None
    assert sphere.surface_id is None
    assert plane.surface_id is None
    assert cylinder.surface_id is None


def test_surface_explicit_id():
    """Test that surface IDs can be explicitly set."""
    # Create surfaces with specific IDs
    sphere = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=1)
    plane = Plane(a=1.0, b=0.0, c=0.0, d=0.0, surface_id=2)
    cylinder = ZCylinder(x0=0.0, y0=0.0, r=1.0, surface_id=3)
    
    # Surface IDs should match what we set
    assert sphere.surface_id == 1
    assert plane.surface_id == 2
    assert cylinder.surface_id == 3


def test_surface_id_setter():
    """Test that surface IDs can be changed after creation."""
    sphere = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0)
    
    # Initially None
    assert sphere.surface_id is None
    
    # Set an ID
    sphere.surface_id = 42
    assert sphere.surface_id == 42
    
    # Change the ID
    sphere.surface_id = 100
    assert sphere.surface_id == 100


def test_geometry_duplicate_surface_ids():
    """Test that geometry validation detects duplicate surface IDs."""
    # Create surfaces with duplicate IDs
    surface1 = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=1)
    surface2 = Plane(a=1.0, b=0.0, c=0.0, d=1.0, surface_id=1)  # Same ID!
    
    # Create regions using these surfaces
    region1 = -surface1  # Inside sphere
    region2 = +surface2  # Above plane
    
    # Create materials and cells
    material1 = Material()
    material1.add_nuclide("H1", 1.0)
    material2 = Material()
    material2.add_nuclide("He4", 1.0)

    cell1 = Cell(name="cell1", region=region1, fill=material1, cell_id=1)
    cell2 = Cell(name="cell2", region=region2, fill=material2, cell_id=2)
    
    # Geometry validation should fail due to duplicate surface IDs
    with pytest.raises(ValueError, match="Duplicate surface_id 1 found"):
        Geometry(cells=[cell1, cell2])


def test_geometry_mixed_surface_ids():
    """Test geometry with mix of surfaces with and without IDs."""
    # Create surfaces - some with IDs, some without
    surface1 = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=1)
    surface2 = Plane(a=1.0, b=0.0, c=0.0, d=1.0)  # No ID
    surface3 = ZCylinder(x0=0.0, y0=0.0, r=0.5, surface_id=3)
    
    # Create regions
    region1 = -surface1 & +surface2  # Inside sphere and above plane
    region2 = -surface3 & -surface2  # Inside cylinder and below plane
    
    # Create materials and cells
    material1 = Material(material_id=10)
    material2 = Material(material_id=20)
    
    cell1 = Cell(name="cell1", region=region1, fill=material1, cell_id=10)
    cell2 = Cell(name="cell2", region=region2, fill=material2, cell_id=20)
    
    # This should work fine - no duplicate surface IDs
    geometry = Geometry(cells=[cell1, cell2])
    assert geometry is not None


def test_geometry_no_surface_ids():
    """Test geometry where no surfaces have IDs."""
    # Create surfaces without IDs
    surface1 = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0)
    surface2 = Plane(a=1.0, b=0.0, c=0.0, d=1.0)
    
    # Create regions
    region1 = -surface1
    region2 = +surface2
    
    # Create materials and cells
    material1 = Material(material_id=1)
    material2 = Material(material_id=2)

    cell1 = Cell(name="cell1", region=region1, fill=material1, cell_id=1)
    cell2 = Cell(name="cell2", region=region2, fill=material2, cell_id=2)

    # This should work fine - no surface ID conflicts
    geometry = Geometry(cells=[cell1, cell2])
    assert geometry is not None


def test_geometry_unique_surface_ids():
    """Test geometry with all surfaces having unique IDs."""
    # Create surfaces with unique IDs
    surface1 = Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=10)
    surface2 = Plane(a=1.0, b=0.0, c=0.0, d=1.0, surface_id=20)
    surface3 = ZCylinder(x0=0.0, y0=0.0, r=0.5, surface_id=30)
    
    # Create regions
    region1 = -surface1 & +surface2
    region2 = -surface3 & -surface2
    
    # Create materials and cells
    material1 = Material(material_id=1)
    material2 = Material(material_id=2)
    
    cell1 = Cell(name="cell1", region=region1, fill=material1, cell_id=1)
    cell2 = Cell(name="cell2", region=region2, fill=material2, cell_id=2)
    
    # This should work fine - all unique surface IDs
    geometry = Geometry(cells=[cell1, cell2])
    assert geometry is not None


def test_surface_types_with_ids():
    """Test that all surface types support optional IDs."""
    surfaces = [
        XPlane(x0=1.0, surface_id=1),
        YPlane(y0=2.0, surface_id=2),
        ZPlane(z0=3.0, surface_id=3),
        Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=4),
        ZCylinder(x0=0.0, y0=0.0, r=0.5, surface_id=5),
        Plane(a=1.0, b=1.0, c=1.0, d=0.0, surface_id=6),
        Cylinder(
            x0=0.0, y0=0.0, z0=0.0,
            axis_x=1.0, axis_y=0.0, axis_z=0.0,
            r=1.0, surface_id=7
        )
    ]
    
    # Check that all surfaces have the expected IDs
    for i, surface in enumerate(surfaces, 1):
        assert surface.surface_id == i


def test_surface_types_without_ids():
    """Test that all surface types work without IDs."""
    surfaces = [
        XPlane(x0=1.0),
        YPlane(y0=2.0),
        ZPlane(z0=3.0),
        Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0),
        ZCylinder(x0=0.0, y0=0.0, r=0.5),
        Plane(a=1.0, b=1.0, c=1.0, d=0.0),
        Cylinder(
            x0=0.0, y0=0.0, z0=0.0,
            axis_x=1.0, axis_y=0.0, axis_z=0.0,
            r=1.0
        )
    ]
    
    # Check that all surfaces have None IDs
    for surface in surfaces:
        assert surface.surface_id is None


def test_complex_geometry_surface_validation():
    """Test surface ID validation in a more complex geometry."""
    # Create multiple surfaces with some duplicate IDs
    surfaces = [
        Sphere(x0=0.0, y0=0.0, z0=0.0, r=2.0, surface_id=1),
        ZCylinder(x0=0.0, y0=0.0, r=1.0, surface_id=2),
        ZPlane(z0=1.0, surface_id=3),
        ZPlane(z0=-1.0, surface_id=3),  # Duplicate ID!
        XPlane(x0=0.0),  # No ID
    ]
    
    # Create a complex region using all surfaces
    region = (-surfaces[0] & -surfaces[1] & +surfaces[2] & -surfaces[3]) | +surfaces[4]
    
    # Create material and cell
    material = Material(material_id=1)
    cell = Cell(name="test_cell", region=region, fill=material, cell_id=1)

    # Should fail due to duplicate surface ID 3
    with pytest.raises(ValueError, match="Duplicate surface_id 3 found"):
        Geometry(cells=[cell])

