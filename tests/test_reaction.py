def test_reaction_instantiation():
    from materials_for_mc import Reaction
    # Provide all required arguments for Reaction
    from materials_for_mc import ReactionProduct
    products = [
        ReactionProduct(
            "n", "emission", 0.0, [], []
        ),
        ReactionProduct(
            "n", "emission", 0.0, [], []
        )
    ]
    r = Reaction(["n"], products, 1.23, [1.23], 0, [0], 18, [0.0, 1.0, 2.0])
    assert r.reactants == ["n"]
    assert [p.particle for p in r.products] == ["n", "n"]
    assert r.energy == 1.23
