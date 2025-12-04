import yaml as mc

def test_python_source_construction():
    # Test new OpenMC-compatible API
    src = mc.IndependentSource()
    src.space = [1.0, 2.0, 3.0]
    src.angle = mc.stats.Monodirectional([0.0, 0.0, 1.0])
    src.energy = 2e6
    
    particle = src.sample()
    assert particle.position == [1.0, 2.0, 3.0]
    assert particle.direction == [0.0, 0.0, 1.0]
    assert particle.energy == 2e6
