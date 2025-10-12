import materials_for_mc as mc


def test_python_settings_construction():
    # Use new OpenMC-compatible API
    src = mc.IndependentSource()
    src.space = [0.0, 0.0, 0.0]
    src.angle = mc.stats.Monodirectional([0.0, 1.0, 0.0])
    src.energy = 1e5
    
    settings = mc.Settings(100, 10, src)
    assert settings.particles == 100
    assert settings.batches == 10
    assert settings.source.space == [0.0, 0.0, 0.0]
    assert settings.source.energy == 14.06e6

def test_python_settings_source_assignment():
    # Use new OpenMC-compatible API
    src = mc.IndependentSource(
        space=[1.0, 1.0, 1.0],
        angle=mc.stats.Monodirectional([0.0, 0.0, 1.0]),
        energy=12e6
    )
    
    settings = mc.Settings(50, 5, src)
    settings.source = src
    assert settings.source.space == [1.0, 1.0, 1.0]
    assert settings.source.energy == 12e6

def test_python_settings_source_assignment():
    # Use new OpenMC-compatible API
    src = mc.IndependentSource(
        space=[1.0, 1.0, 1.0],
        angle=mc.stats.Isotropic(),
        energy=1e6
    )
    
    settings = mc.Settings(50, 5, src)
    settings.source = src
    assert settings.source.space == [1.0, 1.0, 1.0]
    assert settings.source.energy == 1e6
