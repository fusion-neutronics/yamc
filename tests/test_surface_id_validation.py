"""
Test surface ID validation functionality.
Tests optional surface IDs and duplicate surface ID detection in geometry validation.
"""

import pytest
import materials_for_mc as mc


def test_surface_optional_id():
    """Test that surface IDs are optional and default to None."""
    # Create surface without specifying ID
    surface = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0)
    assert surface.surface_id is None
    
    # Create surface with explicit None
    surface2 = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=None)
    assert surface2.surface_id is None


def test_surface_explicit_id():
    """Test that surface IDs can be set explicitly."""
    surface = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=42)
    assert surface.surface_id == 42
    
    # Test with different surface types
    plane = mc.Plane(a=1.0, b=0.0, c=0.0, d=1.0, surface_id=100)
    assert plane.surface_id == 100


def test_surface_id_setter():
    """Test that surface IDs can be modified after creation."""
    surface = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0)
    assert surface.surface_id is None
    
    # Set ID
    surface.surface_id = 123
    assert surface.surface_id == 123
    
    # Set back to None
    surface.surface_id = None
    assert surface.surface_id is None


def test_geometry_duplicate_surface_ids():
    """Test that duplicate surface IDs are detected in geometry validation."""
    # Create surfaces with same ID
    surface1 = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=10)
    surface2 = mc.Sphere(x0=2.0, y0=0.0, z0=0.0, r=1.0, surface_id=10)  # Same ID
    
    # Use different materials to avoid duplicate material ID error
    material1 = mc.Material(material_id=1)
    material2 = mc.Material(material_id=2)
    
    # Create cells using these surfaces
    cell1 = mc.Cell(name='cell1', region=-surface1, fill=material1, cell_id=1)
    cell2 = mc.Cell(name='cell2', region=-surface2, fill=material2, cell_id=2)
    
    # Should raise error due to duplicate surface IDs
    with pytest.raises(ValueError, match="Duplicate surface_id 10 found"):
        mc.Geometry(cells=[cell1, cell2])


def test_geometry_mixed_surface_ids():
    """Test geometry with mix of surfaces with and without IDs."""
    # Create surfaces - mix of with/without IDs
    surface1 = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=10)  # With ID
    surface2 = mc.Sphere(x0=2.0, y0=0.0, z0=0.0, r=1.0)  # Without ID
    surface3 = mc.Sphere(x0=4.0, y0=0.0, z0=0.0, r=1.0, surface_id=20)  # With ID
    
    # Use different materials to avoid duplicate material ID error
    material1 = mc.Material(material_id=1)
    material2 = mc.Material(material_id=2)
    material3 = mc.Material(material_id=3)
    
    # Create cells
    cell1 = mc.Cell(name='cell1', region=-surface1, fill=material1, cell_id=1)
    cell2 = mc.Cell(name='cell2', region=-surface2, fill=material2, cell_id=2)
    cell3 = mc.Cell(name='cell3', region=-surface3, fill=material3, cell_id=3)
    
    # Should work fine - no duplicates
    geometry = mc.Geometry(cells=[cell1, cell2, cell3])
    assert len(geometry.cells) == 3


def test_geometry_no_surface_ids():
    """Test geometry with all surfaces having None IDs."""
    # Create surfaces without IDs
    surface1 = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0)
    surface2 = mc.Sphere(x0=2.0, y0=0.0, z0=0.0, r=1.0)
    
    # Use different materials to avoid duplicate material ID error
    material1 = mc.Material(material_id=1)
    material2 = mc.Material(material_id=2)
    
    # Create cells
    cell1 = mc.Cell(name='cell1', region=-surface1, fill=material1, cell_id=1)
    cell2 = mc.Cell(name='cell2', region=-surface2, fill=material2, cell_id=2)
    
    # Should work fine - None IDs are allowed
    geometry = mc.Geometry(cells=[cell1, cell2])
    assert len(geometry.cells) == 2


def test_geometry_unique_surface_ids():
    """Test geometry with all unique surface IDs."""
    # Create surfaces with unique IDs
    surface1 = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=1)
    surface2 = mc.Sphere(x0=2.0, y0=0.0, z0=0.0, r=1.0, surface_id=2)
    surface3 = mc.Sphere(x0=4.0, y0=0.0, z0=0.0, r=1.0, surface_id=3)
    
    # Use different materials to avoid duplicate material ID error
    material1 = mc.Material(material_id=100)
    material2 = mc.Material(material_id=200)
    material3 = mc.Material(material_id=300)
    
    # Create cells
    cell1 = mc.Cell(name='cell1', region=-surface1, fill=material1, cell_id=1)
    cell2 = mc.Cell(name='cell2', region=-surface2, fill=material2, cell_id=2)
    cell3 = mc.Cell(name='cell3', region=-surface3, fill=material3, cell_id=3)
    
    # Should work fine
    geometry = mc.Geometry(cells=[cell1, cell2, cell3])
    assert len(geometry.cells) == 3


def test_surface_types_with_ids():
    """Test that all surface types support surface_id parameter."""
    # Test various surface types with IDs
    sphere = mc.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=1)
    plane = mc.Plane(a=1, b=0, c=0, d=1, surface_id=2)
    xplane = mc.XPlane(x0=1, surface_id=3)
    yplane = mc.YPlane(y0=2, surface_id=4)
    zplane = mc.ZPlane(z0=3, surface_id=5)
    zcylinder = mc.ZCylinder(x0=0, y0=0, r=1, surface_id=6)
    cylinder = mc.Cylinder(x0=0, y0=0, z0=0, axis_x=0, axis_y=0, axis_z=1, r=1, surface_id=7)
    
    assert sphere.surface_id == 1
    assert plane.surface_id == 2
    assert xplane.surface_id == 3
    assert yplane.surface_id == 4
    assert zplane.surface_id == 5
    assert zcylinder.surface_id == 6
    assert cylinder.surface_id == 7


def test_surface_types_without_ids():
    """Test that all surface types default to None for surface_id."""
    # Test various surface types without IDs
    sphere = mc.Sphere(x0=0, y0=0, z0=0, r=1)
    plane = mc.Plane(a=1, b=0, c=0, d=1)
    xplane = mc.XPlane(x0=1)
    yplane = mc.YPlane(y0=2)
    zplane = mc.ZPlane(z0=3)
    zcylinder = mc.ZCylinder(x0=0, y0=0, r=1)
    cylinder = mc.Cylinder(x0=0, y0=0, z0=0, axis_x=0, axis_y=0, axis_z=1, r=1)
    
    assert sphere.surface_id is None
    assert plane.surface_id is None
    assert xplane.surface_id is None
    assert yplane.surface_id is None
    assert zplane.surface_id is None
    assert zcylinder.surface_id is None
    assert cylinder.surface_id is None


def test_complex_geometry_surface_validation():
    """Test surface validation in a complex geometry with shared surfaces."""
    # Create surfaces with unique IDs
    surface1 = mc.Sphere(x0=0.0, y0=0.0, z0=0.0, r=1.0, surface_id=10)
    surface2 = mc.Plane(a=1.0, b=0.0, c=0.0, d=1.0, surface_id=20)
    surface3 = mc.ZCylinder(x0=0.0, y0=0.0, r=0.5, surface_id=30)
    
    material1 = mc.Material(material_id=1)
    material2 = mc.Material(material_id=2)
    
    # Create regions that share surfaces (this tests Arc deduplication)
    region1 = -surface1 & +surface2   # Uses surface1(10) and surface2(20)
    region2 = -surface3 & -surface2   # Uses surface3(30) and surface2(20) - SHARED!
    
    cell1 = mc.Cell(name='cell1', region=region1, fill=material1, cell_id=1)
    cell2 = mc.Cell(name='cell2', region=region2, fill=material2, cell_id=2)
    
    # Should work fine - shared surfaces are correctly deduplicated
    geometry = mc.Geometry(cells=[cell1, cell2])
    assert len(geometry.cells) == 2


if __name__ == "__main__":
    pytest.main([__file__])
