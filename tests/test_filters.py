import pytest
import materials_for_mc as mmc


def test_cell_filter_creation():
    """Test that CellFilter can be created from a cell."""
    # Create a simple sphere surface
    sphere = mmc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=2.0, boundary_type='vacuum')
    region = -sphere
    
    # Create a material
    material = mmc.Material()
    material.set_density("g/cc", 1.0)
    material.add_nuclide("Li6", 1.0)
    
    # Create a cell
    cell = mmc.Cell(region, 42, "test_cell", material)
    
    # Create CellFilter - this is the main functionality we need to test
    cell_filter = mmc.CellFilter(cell)
    
    # Verify the filter was created successfully
    assert cell_filter is not None
    assert isinstance(cell_filter, mmc.CellFilter)


def test_cell_filter_with_material_id():
    """Test CellFilter creation with a material that has an ID."""
    # Create surfaces and regions
    sphere = mmc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=2.0, boundary_type='vacuum')
    region = -sphere
    
    # Create material with material_id
    material = mmc.Material(name="test_material", material_id=123)
    material.set_density("g/cc", 1.0)
    material.add_nuclide("Li6", 1.0)
    
    # Create cell
    cell = mmc.Cell(region, 42, "test_cell_with_material_id", material)
    
    # Create CellFilter
    cell_filter = mmc.CellFilter(cell)
    assert cell_filter is not None


def test_cell_filter_invalid_input():
    """Test CellFilter behavior with invalid inputs."""
    # Test with None (should raise error during construction)
    with pytest.raises((TypeError, AttributeError)):
        mmc.CellFilter(None)
