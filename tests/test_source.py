import materials_for_mc as mc

def test_python_source_construction():
    src = mc.Source([1.0, 2.0, 3.0], [0.0, 0.0, 1.0], 2e6)
    pos, dir, en = src.sample()
    assert pos == [1.0, 2.0, 3.0]
    assert dir == [0.0, 0.0, 1.0]
    assert en == 2e6
