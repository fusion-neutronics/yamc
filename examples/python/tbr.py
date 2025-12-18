import time


# for code in ['openmc', 'yamc']:
for code in ['yamc']:
    if code == 'openmc':
        import openmc as mc
        # mc.config['cross_sections'] = "/home/jon/nuclear_data/tendl-2021-hdf5/tendl-2021-hdf5/cross_sections.xml"
        mc.config['cross_sections'] = "/home/jon/nuclear_data/cross_sections.xml"
    elif code == 'yamc':
        import yamc as mc

    
    # Create two-cell geometry: inner sphere (Li6) and outer annular region (Be9)
    sphere1 = mc.Sphere(
        surface_id=1,
        x0=0.0,
        y0=0.0,
        z0=0.0,
        r=1.0,
        boundary_type='transmission',
    )
    sphere2 = mc.Sphere(
        surface_id=2,
        x0=0.0,
        y0=0.0,
        z0=0.0,
        r=200.0,
        boundary_type='vacuum',
    )
    region1 = -sphere1
    region2 = +sphere1 & -sphere2

    # Create materials with different absorption characteristics
    material1 = mc.Material()
    material1.material_id = 1  # Set material_id for MaterialFilter testing
    material1.add_nuclide("Li6", 0.07/2)  # Li4SiO4
    material1.add_nuclide("Li7", 0.93/2)
    material1.add_nuclide("Be9", 0.5)
    # material1.add_element("O", 4.0)
    # material1.add_element("Si", 1.0)
    material1.set_density("g/cm3", 2.0)
    if code == 'yamc':
        material1.read_nuclides_from_json({"Be9": "tests/Be9.h5", "Li6": "tests/Li6.h5", "Li7": "tests/Li7.h5"})


    # Create cells
    cell1 = mc.Cell(
        cell_id=1,
        name="inner_sphere",
        region=region1,
    )
    cell2 = mc.Cell(
        cell_id=2,
        name="outer_annular",
        region=region2,
        fill=material1,
    )
    geometry = mc.Geometry([cell1, cell2])

    if code == 'openmc':
        source = mc.IndependentSource(
            space=mc.stats.Point([0,0,0]),
            angle=mc.stats.Isotropic(),
            energy=mc.stats.Discrete([14060000.0], [1.0])
        )
    elif code == 'yamc':
        source = mc.IndependentSource(
            space=[0.0, 0.0, 0.0],
            angle=mc.stats.Isotropic(),
            energy=14.06e6
        )

    settings = mc.Settings(
        particles=50000,
        batches=5,
        source=source,
        seed=1,
    )
    if code == 'openmc':
        settings.run_mode='fixed source'


    # Create tallies with CellFilters
    cell_filter2 = mc.CellFilter(cell2)
    tally1 = mc.Tally()
    tally1.filters = [cell_filter2]
    if code == 'openmc':
        tally1.scores = ['105']  # n,t
    elif code == 'yamc':
        tally1.scores = [105]  # n,t
    tally1.name = "tbr"
    tallies = [tally1]

    model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)

    time.start = time.time()
    if code == 'openmc':
        model.run(apply_tally_results=True)
    elif code == 'yamc':
        model.run()
    print(f"Simulation completed in {time.time() - time.start} seconds.")
    # Tallies are updated in place!
    print(f"{code} TBR (tritium breeding ratio): {tally1.mean}")
