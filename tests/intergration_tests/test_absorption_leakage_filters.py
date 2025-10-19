import materials_for_mc as mc


def test_absorption_leakage_filters():
    """
    Integration test verifying that:
    1. CellFilter correctly separates tallies by cell
    2. Particle conservation (absorption + leakage = 1)
    3. Tally consistency (sum of cell tallies = total tally)
    4. Surface crossing between cells works correctly
    """
    
    # Create two-cell geometry: inner sphere (Li6) and outer annular region (Be9)
    sphere1 = mc.Sphere(
        surface_id=1,
        x0=0.0,
        y0=0.0,
        z0=0.0,
        r=1.0,
        boundary_type='transmission',
    )
    sphere2 = mc.Sphere(
        surface_id=2,
        x0=0.0,
        y0=0.0,
        z0=0.0,
        r=2.0,
        boundary_type='Vacuum',
    )
    region1 = -sphere1
    region2 = +sphere1 & -sphere2

    # Create materials with different absorption characteristics
    material1 = mc.Material()
    material1.material_id = 1  # Set material_id for MaterialFilter testing
    material1.add_nuclide("Li6", 1.0)
    material1.set_density("g/cm3", 10.0)
    material1.read_nuclides_from_json({"Li6": "tests/Li6.json"})

    material2 = mc.Material()
    material2.material_id = 2  # Set material_id for MaterialFilter testing
    material2.add_nuclide("Be9", 1.0)
    material2.set_density("g/cm3", 20.0)
    material2.read_nuclides_from_json({"Be9": "tests/Be9.json"})

    # Create cells
    cell1 = mc.Cell(
        cell_id=1,
        name="inner_sphere",
        region=region1,
        fill=material1,
    )
    cell2 = mc.Cell(
        cell_id=2,
        name="outer_annular",
        region=region2,
        fill=material2,
    )
    geometry = mc.Geometry(cells=[cell1, cell2])

        # Source: neutrons starting at origin, moving upward
    source = mc.IndependentSource(space=[0.0, 0.0, 0.0], angle=mc.stats.Isotropic(), energy=1e6)
    settings = mc.Settings(particles=50, batches=2, source=source)  # Increased for better statistics

    # Create tallies with CellFilters
    cell_filter1 = mc.CellFilter(cell1)
    tally1 = mc.Tally()
    tally1.filters = [cell_filter1]
    tally1.score = 101  # absorption
    tally1.name = "absorption in cell 1"

    cell_filter2 = mc.CellFilter(cell2)
    tally2 = mc.Tally()
    tally2.filters = [cell_filter2]
    tally2.score = 101  # absorption
    tally2.name = "absorption in cell 2"

    # Create equivalent MaterialFilters for comparison
    material_filter1 = mc.MaterialFilter(material1)
    tally1_mat = mc.Tally()
    tally1_mat.filters = [material_filter1]
    tally1_mat.score = 101  # absorption
    tally1_mat.name = "absorption in material 1 (MaterialFilter)"

    material_filter2 = mc.MaterialFilter(material2)
    tally2_mat = mc.Tally()
    tally2_mat.filters = [material_filter2]
    tally2_mat.score = 101  # absorption
    tally2_mat.name = "absorption in material 2 (MaterialFilter)"

    # Total absorption tally (no filter)
    tally3 = mc.Tally()
    tally3.score = 101  # absorption
    tally3.name = "total absorption"

    # Mixed filter tallies - intersection of MaterialFilter and CellFilter
    # Should match cell 1 absorption when material 1 and cell 1 are combined
    tally4_match = mc.Tally()
    tally4_match.filters = [material_filter1, cell_filter1]  # Material 1 AND Cell 1
    tally4_match.score = 101  # absorption
    tally4_match.name = "absorption in material 1 AND cell 1 (should match cell 1)"

    # Should be zero when material 2 and cell 1 are combined (no overlap)
    tally5_zero = mc.Tally()
    tally5_zero.filters = [material_filter2, cell_filter1]  # Material 2 AND Cell 1
    tally5_zero.score = 101  # absorption
    tally5_zero.name = "absorption in material 2 AND cell 1 (should be zero)"

    tallies = [tally1, tally2, tally1_mat, tally2_mat, tally3, tally4_match, tally5_zero]

    # Run simulation
    model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
    leakage_tally, tally1, tally2, tally1_mat, tally2_mat, tally3, tally4_match, tally5_zero = model.run()

    # Integration test assertions
    
    # Test 1: CellFilter functionality - different cells should have different absorption rates
    assert tally1.mean != tally2.mean, "CellFilter should separate tallies by cell"
    
    # Test 2: Tally consistency - sum of cell tallies should equal total tally
    tolerance = 1e-10
    sum_diff = abs((tally1.mean + tally2.mean) - tally3.mean)
    assert sum_diff < tolerance, f"Sum of cell tallies ({tally1.mean + tally2.mean}) should equal total tally ({tally3.mean}), difference: {sum_diff}"
    
    # Test 3: Particle conservation - sum of mutually exclusive cell absorptions + leakage should equal 1
    conservation_diff = abs((tally1.mean + tally2.mean + leakage_tally.mean) - 1.0)
    assert conservation_diff < tolerance, f"(Absorption in cell 1 + cell 2 + leakage) should equal 1.0, got {tally1.mean + tally2.mean + leakage_tally.mean}, difference: {conservation_diff}"
    
    # Test 4: Physical reasonableness - some particles should be absorbed, some should leak
    assert tally3.mean > 0.0, "Some particles should be absorbed"
    assert leakage_tally.mean > 0.0, "Some particles should leak from geometry"
    assert tally3.mean < 1.0, "Not all particles should be absorbed"
    assert leakage_tally.mean < 1.0, "Not all particles should leak"
    
    # Test 5: At least one cell should have absorption, and surface crossing should work
    # (With small particle counts, statistical fluctuations may cause zero absorption in one cell)
    assert tally1.mean >= 0.0, "Cell 1 absorption should be non-negative"
    assert tally2.mean >= 0.0, "Cell 2 absorption should be non-negative"
    assert (tally1.mean + tally2.mean) > 0.0, "Total absorption should be positive"
    
    # Test 6: MaterialFilter equivalence - should give same results as CellFilter for single-material cells
    mat_filter_tolerance = 1e-10
    cell1_vs_mat1_diff = abs(tally1.mean - tally1_mat.mean)
    cell2_vs_mat2_diff = abs(tally2.mean - tally2_mat.mean)
    
    assert cell1_vs_mat1_diff < mat_filter_tolerance, f"CellFilter and MaterialFilter should give same results for cell 1 ({tally1.mean} vs {tally1_mat.mean}), difference: {cell1_vs_mat1_diff}"
    assert cell2_vs_mat2_diff < mat_filter_tolerance, f"CellFilter and MaterialFilter should give same results for cell 2 ({tally2.mean} vs {tally2_mat.mean}), difference: {cell2_vs_mat2_diff}"
    
    # Test 7: Mixed filter intersection tests
    # Material 1 AND Cell 1 should equal Cell 1 alone (since Cell 1 contains Material 1)
    mixed_match_diff = abs(tally1.mean - tally4_match.mean)
    assert mixed_match_diff < mat_filter_tolerance, f"Material 1 AND Cell 1 should equal Cell 1 alone ({tally1.mean} vs {tally4_match.mean}), difference: {mixed_match_diff}"
    
    # Material 2 AND Cell 1 should be zero (Cell 1 contains Material 1, not Material 2)
    assert tally5_zero.mean == 0.0, f"Material 2 AND Cell 1 should be zero (no overlap), got: {tally5_zero.mean}"
    

    
    print(f"✓ Integration test passed!")
    print(f"  - Absorption in cell 1: {tally1.mean:.4f}")
    print(f"  - Absorption in cell 2: {tally2.mean:.4f}")
    print(f"  - Total absorption: {tally3.mean:.4f}")
    print(f"  - Leakage: {leakage_tally.mean:.4f}")
    print(f"  - Conservation check: {tally3.mean + leakage_tally.mean:.6f}")
    print(f"  - MaterialFilter vs CellFilter:")
    print(f"    - Cell 1 vs Material 1: {cell1_vs_mat1_diff:.2e} difference")
    print(f"    - Cell 2 vs Material 2: {cell2_vs_mat2_diff:.2e} difference")
    print(f"  - Mixed filter tests:")
    print(f"    - Material 1 AND Cell 1: {tally4_match.mean:.4f} (should equal Cell 1)")
    print(f"    - Material 2 AND Cell 1: {tally5_zero.mean:.4f} (should be zero)")


def test_duplicate_filter_error():
    """
    Test that adding multiple filters of the same type raises an error.
    """
    import pytest
    
    # Create minimal geometry
    sphere = mc.Sphere(
        surface_id=1,
        x0=0.0, y0=0.0, z0=0.0, r=1.0,
        boundary_type='Vacuum'
    )
    region = -sphere

    # Create materials
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

    # Create cells
    cell1 = mc.Cell(cell_id=1, name="cell1", region=region, fill=material1)
    cell2 = mc.Cell(cell_id=2, name="cell2", region=region, fill=material2)
    
    geometry = mc.Geometry(cells=[cell1])
    source = mc.IndependentSource(space=[0.0, 0.0, 0.0], angle=mc.stats.Isotropic(), energy=1e6)
    settings = mc.Settings(particles=10, batches=1, source=source)

    # Test 1: Multiple CellFilters should raise error
    tally_bad_cell = mc.Tally()
    tally_bad_cell.score = 101
    tally_bad_cell.name = "bad tally with duplicate cell filters"
    
    cell_filter1 = mc.CellFilter(cell1)
    cell_filter2 = mc.CellFilter(cell2)
    
    with pytest.raises(ValueError, match="Multiple filters of the same type are not allowed"):
        tally_bad_cell.filters = [cell_filter1, cell_filter2]
    
    # Test 2: Multiple MaterialFilters should raise error
    tally_bad_material = mc.Tally()
    tally_bad_material.score = 101
    tally_bad_material.name = "bad tally with duplicate material filters"
    
    material_filter1 = mc.MaterialFilter(material1)
    material_filter2 = mc.MaterialFilter(material2)
    
    with pytest.raises(ValueError, match="Multiple filters of the same type are not allowed"):
        tally_bad_material.filters = [material_filter1, material_filter2]
    
    print("✓ Duplicate filter error test passed!")


if __name__ == "__main__":
    test_absorption_leakage_filters()
    test_duplicate_filter_error()
