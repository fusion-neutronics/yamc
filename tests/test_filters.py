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


def test_material_filter_creation():
    """Test that MaterialFilter can be created from a material with ID."""
    # Create a material with ID
    material = mmc.Material(name="test_material", material_id=123)
    material.set_density("g/cc", 1.0)
    material.add_nuclide("Li6", 1.0)
    
    # Create MaterialFilter - this is the main functionality we need to test
    material_filter = mmc.MaterialFilter(material)
    
    # Verify the filter was created successfully
    assert material_filter is not None
    assert isinstance(material_filter, mmc.MaterialFilter)


def test_material_filter_creation_no_id_fails():
    """Test MaterialFilter creation fails with a material that has no ID."""
    # Create material without ID
    material = mmc.Material(name="test_material")
    material.set_density("g/cc", 1.0)
    material.add_nuclide("Li6", 1.0)
    
    # Create MaterialFilter should fail for material with no ID
    with pytest.raises(ValueError, match="Cannot create MaterialFilter for material with no ID"):
        mmc.MaterialFilter(material)


def test_material_filter_representation():
    """Test MaterialFilter string representations."""
    # Test with material that has ID
    material_with_id = mmc.Material(name="test_material", material_id=123)
    filter_with_id = mmc.MaterialFilter(material_with_id)
    
    repr_str = repr(filter_with_id)
    assert "MaterialFilter" in repr_str
    assert "123" in repr_str
    
    str_str = str(filter_with_id)
    assert "MaterialFilter" in str_str


def test_tally_with_material_filter():
    """Test that Tally can use MaterialFilter."""
    material = mmc.Material(name="test_material", material_id=123)
    material.set_density("g/cc", 1.0)
    material.add_nuclide("Li6", 1.0)
    
    # Create tally with MaterialFilter
    tally = mmc.Tally()
    material_filter = mmc.MaterialFilter(material)
    tally.filters = [material_filter]
    tally.scores = [101]  # absorption
    
    # Verify tally has the filter
    assert len(tally.filters) == 1
    # The filter should be returned as the appropriate Python object
    retrieved_filter = tally.filters[0]
    assert isinstance(retrieved_filter, mmc.MaterialFilter)


def test_tally_with_mixed_filters():
    """Test that Tally can use both CellFilter and MaterialFilter together."""
    # Create cell with CellFilter
    sphere = mmc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=2.0, boundary_type='vacuum')
    region = -sphere
    material1 = mmc.Material(material_id=1)
    material1.set_density("g/cc", 1.0)
    material1.add_nuclide("Li6", 1.0)
    cell = mmc.Cell(region, 42, "test_cell", material1)
    
    # Create material with MaterialFilter
    material2 = mmc.Material(name="test_material", material_id=123)
    material2.set_density("g/cc", 1.0)
    material2.add_nuclide("Be9", 1.0)
    
    # Create tally with mixed filters
    tally = mmc.Tally()
    cell_filter = mmc.CellFilter(cell)
    material_filter = mmc.MaterialFilter(material2)
    tally.filters = [cell_filter, material_filter]
    tally.scores = [101]  # absorption
    
    # Verify tally has both filters
    assert len(tally.filters) == 2
    retrieved_filters = tally.filters
    
    # Check that we got the right types back
    filter_types = [type(f).__name__ for f in retrieved_filters]
    assert "CellFilter" in filter_types
    assert "MaterialFilter" in filter_types
