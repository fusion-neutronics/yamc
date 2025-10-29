import materials_for_mc as mc


def test_reproducibility_with_same_seed():
    """Test that simulations with the same seed produce identical results"""
    # Create simple geometry
    sphere = mc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=2.0, boundary_type='vacuum')
    region = -sphere

    # Create material
    material = mc.Material()
    material.material_id = 1
    material.add_nuclide("Li6", 1.0)
    material.set_density("g/cm3", 10.0)
    material.read_nuclides_from_json({"Li6": "tests/Li6.json"})

    # Create cell and geometry
    cell = mc.Cell(cell_id=1, name="test_cell", region=region, fill=material)
    geometry = mc.Geometry(cells=[cell])

    # Source and settings with seed
    source = mc.IndependentSource(space=[0.0, 0.0, 0.0], angle=mc.stats.Isotropic(), energy=1e6)
    settings = mc.Settings(particles=100, batches=10, source=source, seed=42)

    # Create separate tallies for each run
    tally1 = mc.Tally()
    tally1.score = 101  # absorption
    tally1.name = "test_absorption_1"

    tally2 = mc.Tally()
    tally2.score = 101  # absorption
    tally2.name = "test_absorption_2"

    tally3 = mc.Tally()
    tally3.score = 101  # absorption
    tally3.name = "test_absorption_3"

    # Run simulation 1
    model1 = mc.Model(geometry=geometry, settings=settings, tallies=[tally1])
    model1.run()

    # Run simulation 2 with same seed
    model2 = mc.Model(geometry=geometry, settings=settings, tallies=[tally2])
    model2.run()

    # Run simulation 3 with same seed
    model3 = mc.Model(geometry=geometry, settings=settings, tallies=[tally3])
    model3.run()

    # Verify all three runs produced identical results
    # Tallies are updated in place - access results directly

    # Check absorption (no more leakage tally)
    assert tally1.mean == tally2.mean, f"Absorption should be identical with same seed: {tally1.mean} vs {tally2.mean}"
    assert tally1.mean == tally3.mean, f"Absorption should be identical with same seed: {tally1.mean} vs {tally3.mean}"
    assert tally1.batch_data == tally2.batch_data, "Absorption batch data should be identical"
    assert tally1.batch_data == tally3.batch_data, "Absorption batch data should be identical"

    print(f"✓ Python reproducibility test passed!")
    print(f"  Run 1 - Absorption: {tally1.mean:.6f}")
    print(f"  Run 2 - Absorption: {tally2.mean:.6f}")
    print(f"  Run 3 - Absorption: {tally3.mean:.6f}")


def test_different_seeds_produce_different_results():
    """Test that simulations with different seeds produce different results"""
    # Create simple geometry
    sphere = mc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=2.0, boundary_type='vacuum')
    region = -sphere

    # Create material
    material = mc.Material()
    material.material_id = 1
    material.add_nuclide("Li6", 1.0)
    material.set_density("g/cm3", 10.0)
    material.read_nuclides_from_json({"Li6": "tests/Li6.json"})

    # Create cell and geometry
    cell = mc.Cell(cell_id=1, name="test_cell", region=region, fill=material)
    geometry = mc.Geometry(cells=[cell])

    # Source
    source = mc.IndependentSource(space=[0.0, 0.0, 0.0], angle=mc.stats.Isotropic(), energy=1e6)

    # Settings with different seeds
    settings1 = mc.Settings(particles=100, batches=10, source=source, seed=42)
    settings2 = mc.Settings(particles=100, batches=10, source=source, seed=123)

    # Create separate tallies for each run
    tally1 = mc.Tally()
    tally1.score = 101  # absorption
    tally1.name = "test_absorption_1"

    tally2 = mc.Tally()
    tally2.score = 101  # absorption
    tally2.name = "test_absorption_2"

    # Run simulation with seed 42
    model1 = mc.Model(geometry=geometry, settings=settings1, tallies=[tally1])
    model1.run()

    # Run simulation with seed 123
    model2 = mc.Model(geometry=geometry, settings=settings2, tallies=[tally2])
    model2.run()

    # Tallies are updated in place - access results directly
    # Verify different seeds produce different results (with high probability)
    different_absorption = tally1.mean != tally2.mean
    different_batch_data = tally1.batch_data != tally2.batch_data

    assert different_absorption or different_batch_data, \
        f"Different seeds should produce different results (absorption: {tally1.mean} vs {tally2.mean})"

    print(f"✓ Different seeds test passed!")
    print(f"  Seed 42  - Absorption: {tally1.mean:.6f}")
    print(f"  Seed 123 - Absorption: {tally2.mean:.6f}")


if __name__ == "__main__":
    test_reproducibility_with_same_seed()
    test_different_seeds_produce_different_results()
