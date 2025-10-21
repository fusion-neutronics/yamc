import materials_for_mc as mc


def test_reproducibility_with_same_seed():
    """Test that simulations with the same seed produce identical results"""
    # Create simple geometry
    sphere = mc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=2.0, boundary_type='Vacuum')
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

    # Create tally
    tally = mc.Tally()
    tally.score = 101  # absorption
    tally.name = "test_absorption"

    # Run simulation 1
    model1 = mc.Model(geometry=geometry, settings=settings, tallies=[tally])
    results1 = model1.run()

    # Run simulation 2 with same seed
    model2 = mc.Model(geometry=geometry, settings=settings, tallies=[tally])
    results2 = model2.run()

    # Run simulation 3 with same seed
    model3 = mc.Model(geometry=geometry, settings=settings, tallies=[tally])
    results3 = model3.run()

    # Verify all three runs produced identical results
    leakage1, absorption1 = results1
    leakage2, absorption2 = results2
    leakage3, absorption3 = results3

    # Check leakage
    assert leakage1.mean == leakage2.mean, f"Leakage should be identical with same seed: {leakage1.mean} vs {leakage2.mean}"
    assert leakage1.mean == leakage3.mean, f"Leakage should be identical with same seed: {leakage1.mean} vs {leakage3.mean}"
    assert leakage1.batch_data == leakage2.batch_data, "Leakage batch data should be identical"
    assert leakage1.batch_data == leakage3.batch_data, "Leakage batch data should be identical"

    # Check absorption
    assert absorption1.mean == absorption2.mean, f"Absorption should be identical with same seed: {absorption1.mean} vs {absorption2.mean}"
    assert absorption1.mean == absorption3.mean, f"Absorption should be identical with same seed: {absorption1.mean} vs {absorption3.mean}"
    assert absorption1.batch_data == absorption2.batch_data, "Absorption batch data should be identical"
    assert absorption1.batch_data == absorption3.batch_data, "Absorption batch data should be identical"

    print(f"✓ Python reproducibility test passed!")
    print(f"  Run 1 - Leakage: {leakage1.mean:.6f}, Absorption: {absorption1.mean:.6f}")
    print(f"  Run 2 - Leakage: {leakage2.mean:.6f}, Absorption: {absorption2.mean:.6f}")
    print(f"  Run 3 - Leakage: {leakage3.mean:.6f}, Absorption: {absorption3.mean:.6f}")


def test_different_seeds_produce_different_results():
    """Test that simulations with different seeds produce different results"""
    # Create simple geometry
    sphere = mc.Sphere(surface_id=1, x0=0.0, y0=0.0, z0=0.0, r=2.0, boundary_type='Vacuum')
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

    # Create tally
    tally = mc.Tally()
    tally.score = 101  # absorption
    tally.name = "test_absorption"

    # Run simulation with seed 42
    model1 = mc.Model(geometry=geometry, settings=settings1, tallies=[tally])
    results1 = model1.run()

    # Run simulation with seed 123
    model2 = mc.Model(geometry=geometry, settings=settings2, tallies=[tally])
    results2 = model2.run()

    leakage1, absorption1 = results1
    leakage2, absorption2 = results2

    # Verify different seeds produce different results (with high probability)
    different_leakage = leakage1.mean != leakage2.mean
    different_absorption = absorption1.mean != absorption2.mean
    different_batch_data = absorption1.batch_data != absorption2.batch_data

    assert different_leakage or different_absorption or different_batch_data, \
        f"Different seeds should produce different results (leakage: {leakage1.mean} vs {leakage2.mean}, absorption: {absorption1.mean} vs {absorption2.mean})"

    print(f"✓ Different seeds test passed!")
    print(f"  Seed 42  - Leakage: {leakage1.mean:.6f}, Absorption: {absorption1.mean:.6f}")
    print(f"  Seed 123 - Leakage: {leakage2.mean:.6f}, Absorption: {absorption2.mean:.6f}")


if __name__ == "__main__":
    test_reproducibility_with_same_seed()
    test_different_seeds_produce_different_results()
