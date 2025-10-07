import materials_for_mc as csg
import pytest

def test_geometry_find_cell():
    surf = csg.Sphere(x0=0, y0=0, z0=0, r=2, surface_id=1)
    region = -surf
    cell = csg.Cell(cell_id=1, region=region)
    geometry = csg.Geometry([cell])
    assert geometry.find_cell(0, 0, 0) is not None
    assert geometry.find_cell(5, 0, 0) is None


def test_geometry_duplicate_id_validation():
    """Test geometry validation for duplicate cell and material IDs."""
    
    # Create surfaces for our cells
    surf1 = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=1)
    surf2 = csg.Sphere(x0=2, y0=0, z0=0, r=1, surface_id=2)
    region1 = -surf1
    region2 = -surf2
    
    # Step 1: Create two cells with IDs 1 and 2
    cell1 = csg.Cell(region=region1, cell_id=1)
    cell2 = csg.Cell(region=region2, cell_id=2)
    
    # Step 2: Create one material with ID 1
    material1 = csg.Material(material_id=1)
    material1.add_nuclide("Li6", 1.0)
    
    # Step 3: Fill both cells with the same material - this should fail due to duplicate material IDs
    # (each cell gets its own copy of the material, but with the same ID)
    cell1_filled = csg.Cell(region=region1, cell_id=1, fill=material1)
    cell2_filled = csg.Cell(region=region2, cell_id=2, fill=material1)
    
    # Step 4: Creating geometry should fail due to duplicate material IDs
    with pytest.raises(ValueError, match="Duplicate material_id 1 found"):
        csg.Geometry([cell1_filled, cell2_filled])
    
    # Step 5: Test duplicate cell ID detection
    # Create another cell with the same ID as cell1
    cell3_duplicate_id = csg.Cell(region=region1, cell_id=1)  # Same ID as cell1
    
    # Creating geometry with duplicate cell IDs should fail
    with pytest.raises(ValueError, match="Duplicate cell_id 1 found"):
        csg.Geometry([cell1, cell3_duplicate_id])
    
    # Step 6: Test creating valid geometry with different material IDs
    # Create two different materials with different IDs
    material_a = csg.Material(material_id=10)
    material_a.add_nuclide("Li6", 1.0)
    
    material_b = csg.Material(material_id=20)
    material_b.add_nuclide("Li7", 1.0)
    
    # Fill cells with different materials
    cell1_with_mat_a = csg.Cell(region=region1, cell_id=1, fill=material_a)
    cell2_with_mat_b = csg.Cell(region=region2, cell_id=2, fill=material_b)
    
    # This should work fine - different cell IDs and different material IDs
    geometry_valid = csg.Geometry([cell1_with_mat_a, cell2_with_mat_b])
    assert len(geometry_valid.cells) == 2
    
    # Step 7: Test duplicate material ID detection
    # Create another material with the same ID as material_a
    material_c_duplicate_id = csg.Material(material_id=10)  # Same ID as material_a
    material_c_duplicate_id.add_nuclide("Be9", 1.0)
    
    # Create a cell with the duplicate material
    cell3_with_duplicate_mat = csg.Cell(region=region1, cell_id=3, fill=material_c_duplicate_id)
    
    # Creating geometry with duplicate material IDs should fail
    with pytest.raises(ValueError, match="Duplicate material_id 10 found"):
        csg.Geometry([cell1_with_mat_a, cell3_with_duplicate_mat])
    
    # Step 8: Change cell ID and test that material validation still works
    # Update cell ID to be unique but keep duplicate material ID
    cell3_updated_id = csg.Cell(region=region2, cell_id=30, fill=material_c_duplicate_id)
    
    # Should still fail due to duplicate material IDs
    with pytest.raises(ValueError, match="Duplicate material_id 10 found"):
        csg.Geometry([cell1_with_mat_a, cell3_updated_id])
    
    # Step 9: Fix material ID and verify it works
    material_c_fixed = csg.Material(material_id=30)  # Different ID
    material_c_fixed.add_nuclide("Be9", 1.0)
    
    cell3_fixed = csg.Cell(region=region2, cell_id=30, fill=material_c_fixed)
    
    # Now it should work
    geometry_final = csg.Geometry([cell1_with_mat_a, cell3_fixed])
    assert len(geometry_final.cells) == 2
    
    # Step 10: Test cells without materials are allowed (no material ID conflicts)
    cell_no_material1 = csg.Cell(region=region1, cell_id=100)
    cell_no_material2 = csg.Cell(region=region2, cell_id=200)
    
    geometry_no_materials = csg.Geometry([cell_no_material1, cell_no_material2])
    assert len(geometry_no_materials.cells) == 2


def test_surface_id_validation():
    """Test surface ID validation in geometry creation."""
    
    # Step 1: Test duplicate surface ID detection
    # Create two surfaces with the same ID
    surf1 = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=10)
    surf2 = csg.Sphere(x0=2, y0=0, z0=0, r=1, surface_id=10)  # Same ID - should fail
    
    region1 = -surf1
    region2 = -surf2
    
    cell1 = csg.Cell(region=region1, cell_id=1)
    cell2 = csg.Cell(region=region2, cell_id=2)
    
    # Creating geometry with duplicate surface IDs should fail
    with pytest.raises(ValueError, match="Duplicate surface_id 10 found"):
        csg.Geometry([cell1, cell2])
    
    # Step 2: Test surface without ID validation
    # Create surface without specifying ID (should default to None)
    surf_no_id = csg.Sphere(x0=0, y0=0, z0=0, r=1)  # No surface_id specified
    region_no_id = -surf_no_id
    cell_no_surface_id = csg.Cell(region=region_no_id, cell_id=1)
    
    # Creating geometry with surface without ID should work fine (None IDs are allowed)
    geometry_with_none_id = csg.Geometry([cell_no_surface_id])
    assert len(geometry_with_none_id.cells) == 1
    
    # Step 3: Test valid unique surface IDs
    surf_a = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=1)
    surf_b = csg.Sphere(x0=2, y0=0, z0=0, r=1, surface_id=2)
    surf_c = csg.Sphere(x0=4, y0=0, z0=0, r=1, surface_id=3)
    
    region_a = -surf_a
    region_b = -surf_b
    region_c = -surf_c
    
    cell_a = csg.Cell(region=region_a, cell_id=1)
    cell_b = csg.Cell(region=region_b, cell_id=2)
    cell_c = csg.Cell(region=region_c, cell_id=3)
    
    # This should work fine - all surface IDs are unique
    geometry_valid = csg.Geometry([cell_a, cell_b, cell_c])
    assert len(geometry_valid.cells) == 3
    
    # Step 4: Test different surface types with unique IDs
    sphere_surf = csg.Sphere(x0=0, y0=0, z0=0, r=1, surface_id=100)
    plane_surf = csg.XPlane(x0=5, surface_id=101)
    cylinder_surf = csg.ZCylinder(x0=0, y0=0, r=2, surface_id=102)
    
    sphere_region = -sphere_surf
    plane_region = -plane_surf
    cylinder_region = -cylinder_surf
    
    sphere_cell = csg.Cell(region=sphere_region, cell_id=10)
    plane_cell = csg.Cell(region=plane_region, cell_id=11)
    cylinder_cell = csg.Cell(region=cylinder_region, cell_id=12)
    
    # Different surface types with unique IDs should work
    geometry_mixed = csg.Geometry([sphere_cell, plane_cell, cylinder_cell])
    assert len(geometry_mixed.cells) == 3
    
    # Step 5: Test duplicate surface ID with different surface types
    duplicate_plane = csg.XPlane(x0=10, surface_id=100)  # Same ID as sphere_surf
    duplicate_region = -duplicate_plane
    duplicate_cell = csg.Cell(region=duplicate_region, cell_id=20)
    
    # Should fail due to duplicate surface ID even with different surface types
    with pytest.raises(ValueError, match="Duplicate surface_id 100 found"):
        csg.Geometry([sphere_cell, duplicate_cell])
    
    # Step 6: Test complex regions with multiple surfaces
    # Create a region using intersection of two surfaces with unique IDs
    surf_x = csg.Sphere(x0=0, y0=0, z0=0, r=2, surface_id=200)
    surf_y = csg.XPlane(x0=0, surface_id=201)
    
    # Region inside sphere AND to the right of plane
    complex_region = (-surf_x) & (+surf_y)
    complex_cell = csg.Cell(region=complex_region, cell_id=50)
    
    # Should work fine - both surfaces have unique IDs
    geometry_complex = csg.Geometry([complex_cell])
    assert len(geometry_complex.cells) == 1
    
    # Step 7: Test complex region with duplicate surface IDs
    surf_duplicate = csg.YPlane(y0=0, surface_id=200)  # Same ID as surf_x
    
    # Region using duplicate surface ID
    duplicate_complex_region = (-surf_x) & (+surf_duplicate)
    duplicate_complex_cell = csg.Cell(region=duplicate_complex_region, cell_id=51)
    
    # Should fail due to duplicate surface ID in complex region
    with pytest.raises(ValueError, match="Duplicate surface_id 200 found"):
        csg.Geometry([duplicate_complex_cell])
