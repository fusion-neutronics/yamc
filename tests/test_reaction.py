def test_reaction_instantiation():
    from materials_for_mc import Reaction
    # Provide all required arguments for Reaction
    from materials_for_mc import ReactionProduct
    products = [
        ReactionProduct("neutron", "emission", 0.0),
        ReactionProduct("neutron", "emission", 0.0)
    ]
    r = Reaction(products, [1.23], 0, [0], 18, [0.0, 1.0, 2.0])
    assert [p.particle for p in r.products] == ["neutron", "neutron"]
    assert r.mt_number == 18
    assert r.cross_section == [1.23]
    assert r.energy_grid == [0.0, 1.0, 2.0]
