import materials_for_mc as mc

def test_python_source_construction():
    src = mc.IndependentSource([1.0, 2.0, 3.0], [0.0, 0.0, 1.0], 2e6)
    particle = src.sample()
    assert particle.position == [1.0, 2.0, 3.0]
    assert particle.direction == [0.0, 0.0, 1.0]
    assert particle.energy == 2e6
