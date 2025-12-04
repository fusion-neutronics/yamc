import yamc as mc

def test_model_construction():
    # Create sphere surface with vacuum boundary
    sphere = mc.Sphere(
        surface_id=1,
        x0=0.0,
        y0=0.0,
        z0=0.0,
        r=2.0,
        boundary_type='vacuum',
    )
    region = -sphere
    material = mc.Material()
    material.add_nuclide("Li6", 1.0)
    material.set_density("g/cm3", 0.5)
    material.read_nuclides_from_json({"Li6": "tests/Li6.json"})
    cell = mc.Cell(
        cell_id=1,
        name="sphere_cell",
        region=region,
        fill=material,
    )
    geometry = mc.Geometry(cells=[cell])
    source = mc.IndependentSource(space=[0.0, 0.0, 0.0], angle=mc.stats.Isotropic(), energy=1e6)
    settings = mc.Settings(particles=5, batches=2, source=source)
    model = mc.Model(geometry, settings)
    assert model.settings.particles == 5
    assert model.settings.source.energy == 1e6
    assert len(model.geometry.cells) == 1
    model.run()
