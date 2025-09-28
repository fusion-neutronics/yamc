import materials_for_mc as mc


def test_python_settings_construction():
    src = mc.Source([0.0, 0.0, 0.0], [0.0, 1.0, 0.0], 1e5)
    settings = mc.Settings(100, 10, src)
    assert settings.particles == 100
    assert settings.batches == 10
    assert settings.source.position == [0.0, 0.0, 0.0]

def test_python_settings_source_assignment():
    src = mc.Source([1.0, 1.0, 1.0], [0.0, 0.0, 1.0], 1e6)
    settings = mc.Settings(50, 5, src)
    settings.source = src
    assert settings.source.position == [1.0, 1.0, 1.0]
