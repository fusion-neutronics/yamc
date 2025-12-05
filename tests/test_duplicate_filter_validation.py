import pytest
import yamc as mc


def test_duplicate_filter_validation():
    """
    Test that adding multiple filters of the same type raises a ValueError.
    
    This prevents impossible conditions where a particle would need to be
    in multiple cells or materials simultaneously.
    """
    
    # Create minimal test materials and geometry
    sphere = mc.Sphere(
        surface_id=1,
        x0=0.0, y0=0.0, z0=0.0, r=1.0,
        boundary_type='vacuum'
    )
    region = -sphere

    material1 = mc.Material()
    material1.material_id = 1
    material1.add_nuclide("Li6", 1.0)
    material1.set_density("g/cm3", 10.0)
    material1.read_nuclides_from_json({"Li6": "tests/Li6.json"})

    material2 = mc.Material()
    material2.material_id = 2
    material2.add_nuclide("Be9", 1.0)
    material2.set_density("g/cm3", 20.0)
    material2.read_nuclides_from_json({"Be9": "tests/Be9.json"})

    cell1 = mc.Cell(cell_id=1, name="cell1", region=region, fill=material1)
    cell2 = mc.Cell(cell_id=2, name="cell2", region=region, fill=material2)
    
    # Test 1: Multiple CellFilters should raise ValueError
    tally1 = mc.Tally()
    tally1.scores = [101]  # absorption
    tally1.name = "tally with duplicate cell filters"
    
    cell_filter1 = mc.CellFilter(cell1)
    cell_filter2 = mc.CellFilter(cell2)
    
    with pytest.raises(ValueError) as excinfo:
        tally1.filters = [cell_filter1, cell_filter2]
    
    assert "Multiple filters of the same type are not allowed" in str(excinfo.value)
    assert "CellFilter" in str(excinfo.value)
    
    # Test 2: Multiple MaterialFilters should raise ValueError  
    tally2 = mc.Tally()
    tally2.scores = [101]  # absorption
    tally2.name = "tally with duplicate material filters"
    
    material_filter1 = mc.MaterialFilter(material1)
    material_filter2 = mc.MaterialFilter(material2)
    
    with pytest.raises(ValueError) as excinfo:
        tally2.filters = [material_filter1, material_filter2]
    
    assert "Multiple filters of the same type are not allowed" in str(excinfo.value)
    assert "MaterialFilter" in str(excinfo.value)
    
    # Test 3: Mixed filter types (different types) should be allowed
    tally3 = mc.Tally()
    tally3.scores = [101]  # absorption
    tally3.name = "tally with mixed filter types"
    
    # This should work without raising an exception
    tally3.filters = [cell_filter1, material_filter1]
    assert len(tally3.filters) == 2
    
    # Test 4: Single filter of each type should work
    tally4 = mc.Tally()
    tally4.scores = [101]
    tally4.name = "tally with single cell filter"
    tally4.filters = [cell_filter1]
    assert len(tally4.filters) == 1
    
    tally5 = mc.Tally()
    tally5.scores = [101]
    tally5.name = "tally with single material filter"
    tally5.filters = [material_filter1]
    assert len(tally5.filters) == 1
    

if __name__ == "__main__":
    test_duplicate_filter_validation()
    print("✓ All duplicate filter validation tests passed!")