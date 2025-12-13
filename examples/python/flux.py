import time
import matplotlib.pyplot as plt
import numpy as np

import yamc
# materuals_for_mc.
isotopes = yamc.natural_abundance().keys()
# Plot the energy-binned flux spectrum

# for isotope in isotopes:
for isotope in ['Cr52']:
    fig, ax = plt.subplots(figsize=(10, 6))
    for code in ['openmc', 'yamc_json', 'yamc_h5']:
    # for code in ['yamc']:

        if code == 'openmc':
            import openmc as mc
            mc.config['cross_sections'] = "/home/jon/nuclear_data/tendl-2021-hdf5/tendl-2021-hdf5/cross_sections.xml"
            # mc.config['cross_sections'] = "/home/jon/nuclear_data/cross_sections.xml"
        elif code == 'yamc_json' or code == 'yamc_h5':
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
        material1.add_nuclide(isotope, 01.0)  # Li4SiO4
        material1.set_density("g/cm3", 2.0)
        if code == 'yamc_json':
            dir = "/home/jon/cross_section_data_tendl_2021/tendl_2021/"
            material1.read_nuclides_from_json({isotope: f"{dir}{isotope}.json"})
        if code == 'yamc_h5':
            dir = "/home/jon/nuclear_data/tendl-2021-hdf5/tendl-2021-hdf5/"
            material1.read_nuclides_from_json({isotope: f"{dir}{isotope}.h5"})


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
            settings = mc.Settings(
                particles=50000,
                batches=2,
                source=source,
                seed=1,
                run_mode='fixed source'
            )
        elif code == 'yamc_json' or code == 'yamc_h5':
            source = mc.IndependentSource(
                space=[0.0, 0.0, 0.0],
                angle=mc.stats.Isotropic(),
                energy=14.06e6
            )
            settings = mc.Settings(
                particles=50000,
                batches=2,
                source=source,
                seed=1
            )

        # Create tallies with CellFilters
        cell_filter2 = mc.CellFilter(cell2)
        tally1 = mc.Tally()
        tally1.filters = [cell_filter2]
        tally1.scores = ['flux']
        tally1.name = "flux"

        # Create energy-binned flux tally
        # Energy bins: logarithmically spaced from 0.1 eV to 20 MeV
        energy_bins = np.logspace(np.log10(0.01), np.log10(20e6), 50)  # 49 bins in eV
        energy_filter = mc.EnergyFilter(energy_bins)

        tally2 = mc.Tally()
        tally2.filters = [cell_filter2, energy_filter]
        tally2.scores = ['flux']
        tally2.name = "flux_energy_binned"

        tallies = [tally1, tally2]

        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)

        time.start = time.time()
        if code == 'openmc':
            model.run(apply_tally_results=True)
        elif code == 'yamc_json' or code == 'yamc_h5':
            model.run()

        print(f"Simulation completed in {time.time() - time.start} seconds.")
        # Tallies are updated in place!
        print(f"Total Flux: {tally1.mean}")
        print(f"Energy-binned flux has {len(tally2.mean)} bins")

        if code == 'yamc_json' or code == 'yamc_h5':
            ax.step(energy_bins[:-1], tally2.mean, where='post', label=code)
            bin_centers = np.sqrt(energy_bins[:-1] * energy_bins[1:])
            ax.errorbar(bin_centers, tally2.mean, yerr=tally2.std_dev, fmt='none', capsize=3, color=ax.lines[-1].get_color())
            # assert np.isclose(sum(tally1.mean), tally1.mean[0])
        elif code == 'openmc':
            ax.step(energy_bins[:-1], tally2.mean.squeeze(), where='post', label=code)
            bin_centers = np.sqrt(energy_bins[:-1] * energy_bins[1:])
            ax.errorbar(bin_centers, tally2.mean.squeeze(), yerr=tally2.std_dev.squeeze(), fmt='none', capsize=3, color=ax.lines[-1].get_color())
            # assert np.isclose(sum(tally2.mean.squeeze()), tally2.mean.squeeze())

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Flux (particles/cm²/source particle)')
    ax.set_title(f'Energy-dependent Flux Spectrum {isotope}')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'flux_spectrum_{isotope}.png', dpi=150)
    print(f"Flux spectrum plot saved as 'flux_spectrum_{isotope}.png'")
    # plt.show()
    plt.close()