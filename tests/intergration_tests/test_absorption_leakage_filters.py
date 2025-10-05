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
    material1.add_nuclide("Li6", 1.0)
    material1.set_density("g/cm3", 10.0)
    material1.read_nuclides_from_json({"Li6": "tests/Li6.json"})

    material2 = mc.Material()
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
    source = mc.Source(position=[0.0, 0.0, 0.0], direction=[0.0, 0.0, 1.0], energy=1e6)
    settings = mc.Settings(particles=100, batches=10, source=source)

    # Create tallies with CellFilters
    filter1 = mc.CellFilter(cell1)
    tally1 = mc.Tally()
    tally1.filters = [filter1]
    tally1.score = 101  # absorption
    tally1.name = "absorption in cell 1"

    filter2 = mc.CellFilter(cell2)
    tally2 = mc.Tally()
    tally2.filters = [filter2]
    tally2.score = 101  # absorption
    tally2.name = "absorption in cell 2"

    # Total absorption tally (no filter)
    tally3 = mc.Tally()
    tally3.score = 101  # absorption
    tally3.name = "total absorption"

    tallies = [tally1, tally2, tally3]

    # Run simulation
    model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)
    leakage_tally, tally1, tally2, tally3 = model.run()

    # Integration test assertions
    
    # Test 1: CellFilter functionality - different cells should have different absorption rates
    assert tally1.mean != tally2.mean, "CellFilter should separate tallies by cell"
    
    # Test 2: Tally consistency - sum of cell tallies should equal total tally
    tolerance = 1e-10
    sum_diff = abs((tally1.mean + tally2.mean) - tally3.mean)
    assert sum_diff < tolerance, f"Sum of cell tallies ({tally1.mean + tally2.mean}) should equal total tally ({tally3.mean}), difference: {sum_diff}"
    
    # Test 3: Particle conservation - absorbed + leaked should equal 1
    conservation_diff = abs((tally3.mean + leakage_tally.mean) - 1.0)
    assert conservation_diff < tolerance, f"Absorption + leakage should equal 1.0, got {tally3.mean + leakage_tally.mean}, difference: {conservation_diff}"
    
    # Test 4: Physical reasonableness - some particles should be absorbed, some should leak
    assert tally3.mean > 0.0, "Some particles should be absorbed"
    assert leakage_tally.mean > 0.0, "Some particles should leak from geometry"
    assert tally3.mean < 1.0, "Not all particles should be absorbed"
    assert leakage_tally.mean < 1.0, "Not all particles should leak"
    
    # Test 5: Both cells should contribute to absorption (surface crossing working)
    assert tally1.mean > 0.0, "Cell 1 should have some absorption"
    assert tally2.mean > 0.0, "Cell 2 should have some absorption (proves surface crossing works)"
    
    print(f"✓ Integration test passed!")
    print(f"  - Absorption in cell 1: {tally1.mean:.4f}")
    print(f"  - Absorption in cell 2: {tally2.mean:.4f}")
    print(f"  - Total absorption: {tally3.mean:.4f}")
    print(f"  - Leakage: {leakage_tally.mean:.4f}")
    print(f"  - Conservation check: {tally3.mean + leakage_tally.mean:.6f}")


if __name__ == "__main__":
    test_absorption_leakage_filters()