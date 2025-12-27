import time
import matplotlib.pyplot as plt
import numpy as np
import yamc
import glob
import os

# Path to OpenMC-format HDF5 nuclear data files
# Update this to point to your nuclear data directory


NUCLEAR_DATA_DIR = "/home/jon/nuclear_data/tendl-2021-hdf5/tendl-2021-hdf5/"#fendl-3.2c-hdf5/"
# NUCLEAR_DATA_DIR = "/home/jon/nuclear_data/fendl-3.2c-hdf5/neutron"
# Find all HDF5 files in the neutron data directory and sort them alphabetically

h5_files = sorted(glob.glob(os.path.join(NUCLEAR_DATA_DIR, "*.h5")))

# Prepare to collect average relative differences for all isotopes
isotopes = [os.path.splitext(os.path.basename(f))[0] for f in h5_files]
avg_rel_diff_list = []
isotopes = ['O16']
for isotope in isotopes:
# for isotope in ['S33', 'Er162',  'K39', 'K41', 'Ga71', 'Cd106']:
    # if isotope not in ['O16', 'Cr52']:
    fig, ax = plt.subplots(figsize=(10, 6))

    # Energy bins: logarithmically spaced from 0.01 eV to 20 MeV
    energy_bins = np.logspace(np.log10(0.01), np.log10(20e6), 20)

    # Store results for comparison table
    results = {}

    for code in ['openmc', 'yamc']:
        if code == 'openmc':
            import openmc as mc
            mc.config['cross_sections'] = f"{NUCLEAR_DATA_DIR}cross_sections.xml"
        elif code == 'yamc':
            import yamc as mc

        # Create two-cell geometry: inner sphere and outer annular region
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

        # Create material
        material1 = mc.Material()
        material1.material_id = 1
        material1.add_nuclide(isotope, 1.0)
        material1.set_density("g/cm3", 1.0)
        if code == 'yamc':
            print(f"Reading nuclide data for {isotope}...")
            # Load from OpenMC-format HDF5 file
            material1.read_nuclides_from_h5({isotope: f"{NUCLEAR_DATA_DIR}/{isotope}.h5"})

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
                batches=20,
                source=source,
                seed=1,
                run_mode='fixed source'
            )
        elif code == 'yamc':
            source = mc.IndependentSource(
                space=[0.0, 0.0, 0.0],
                angle=mc.stats.Isotropic(),
                energy=14.06e6
            )
            settings = mc.Settings(
                particles=50000,
                batches=20,
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
        energy_filter = mc.EnergyFilter(energy_bins)

        tally2 = mc.Tally()
        tally2.filters = [cell_filter2, energy_filter]
        tally2.scores = ['flux']
        tally2.name = "flux_energy_binned"

        tallies = [tally1, tally2]

        model = mc.Model(geometry=geometry, settings=settings, tallies=tallies)

        start_time = time.time()
        if code == 'openmc':
            model.run(apply_tally_results=True)
        elif code == 'yamc':
            model.run()

        print(f"{code}: Simulation completed in {time.time() - start_time:.2f} seconds.")
        print(f"{code}: Total Flux: {tally1.mean}")
        print(f"{code}: Energy-binned flux has {len(tally2.mean)} bins")

        # Store results
        if code == 'yamc':
            results['yamc'] = np.array(tally2.mean)
            ax.step(energy_bins[:-1], tally2.mean, where='post', label=code)
            bin_centers = np.sqrt(energy_bins[:-1] * energy_bins[1:])
            ax.errorbar(bin_centers, tally2.mean, yerr=tally2.std_dev, fmt='none', capsize=3, color=ax.lines[-1].get_color())
        elif code == 'openmc':
            results['openmc'] = tally2.mean.squeeze()
            ax.step(energy_bins[:-1], tally2.mean.squeeze(), where='post', label=code)
            bin_centers = np.sqrt(energy_bins[:-1] * energy_bins[1:])
            ax.errorbar(bin_centers, tally2.mean.squeeze(), yerr=tally2.std_dev.squeeze(), fmt='none', capsize=3, color=ax.lines[-1].get_color())

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Energy (eV)')
    ax.set_ylabel('Flux (particles/cm²/source particle)')
    ax.set_title(f'Energy-dependent Flux Spectrum {isotope} {NUCLEAR_DATA_DIR}')
    ax.grid(True, alpha=0.3)
    ax.legend()

    plt.tight_layout()
    plt.savefig(f'flux_spectrum_{isotope}.png', dpi=150)
    print(f"Flux spectrum plot saved as 'flux_spectrum_{isotope}.png'")
    plt.close()

    # Print comparison table of flux values
    print("\n" + "="*110)
    print(f"FLUX COMPARISON TABLE FOR {isotope}")
    print("="*110)
    print(f"{'Bin':<4} {'E_low (eV)':<14} {'E_high (eV)':<14} {'OpenMC':<16} {'YAMC':<16} {'Diff':<16} {'Rel Diff %':<12}")
    print("-"*110)

    openmc_flux = results.get('openmc', np.zeros(len(energy_bins)-1))
    yamc_flux = results.get('yamc', np.zeros(len(energy_bins)-1))


    # Calculate and print bin-by-bin differences, and collect for averaging
    abs_diffs = []
    rel_diffs = []
    for i in range(len(energy_bins) - 1):
        e_low = energy_bins[i]
        e_high = energy_bins[i + 1]
        omc = openmc_flux[i]
        ymc = yamc_flux[i]
        diff = ymc - omc
        if omc != 0:
            rel_diff = diff / omc * 100
        elif ymc != 0:
            rel_diff = float('inf')
        else:
            rel_diff = 0.0

        abs_diffs.append(abs(diff))
        rel_diffs.append(abs(rel_diff) if np.isfinite(rel_diff) else 0.0)

        # Highlight significant differences
        marker = ""
        if abs(rel_diff) > 50 and (omc > 0.01 or ymc > 0.01):
            marker = " ***"
        elif abs(rel_diff) > 20 and (omc > 0.01 or ymc > 0.01):
            marker = " **"
        elif abs(rel_diff) > 10 and (omc > 0.01 or ymc > 0.01):
            marker = " *"

        print(f"{i:<4} {e_low:<14.4e} {e_high:<14.4e} {omc:<16.6e} {ymc:<16.6e} {diff:<+16.6e} {rel_diff:<12.2f}{marker}")


    # Print average difference metrics for this isotope
    avg_abs_diff = np.mean(abs_diffs)
    avg_rel_diff = np.mean(rel_diffs)
    print(f"\nAverage absolute difference per bin: {avg_abs_diff:.4e}")
    print(f"Average relative difference per bin: {avg_rel_diff:.2f}%\n")

    # Store for summary
    avg_rel_diff_list.append((isotope, avg_rel_diff))

    # Save flux spectrum to CSV file (yamc, openmc, energy)
    csv_filename = f'flux_spectrum_{isotope}.txt'
    bin_centers = np.sqrt(energy_bins[:-1] * energy_bins[1:])
    with open(csv_filename, 'w') as f:
        f.write('yamc,openmc,energy\n')
        for i in range(len(bin_centers)):
            f.write(f"{yamc_flux[i]:.8e},{openmc_flux[i]:.8e},{bin_centers[i]:.8e}\n")
    print(f"Flux spectrum data saved as '{csv_filename}' (columns: yamc, openmc, energy)")

    print("-"*110)
    omc_total = np.sum(openmc_flux)
    ymc_total = np.sum(yamc_flux)
    total_diff = ymc_total - omc_total
    total_rel_diff = (total_diff / omc_total * 100) if omc_total != 0 else 0.0
    print(f"{'TOT':<4} {'':<14} {'':<14} {omc_total:<16.6e} {ymc_total:<16.6e} {total_diff:<+16.6e} {total_rel_diff:<12.2f}")
    print("="*110)

    print("\nLegend: * = >10% diff, ** = >20% diff, *** = >50% diff (only for flux > 0.01)")

# After all isotopes, print sorted summary of avg_rel_diff
if avg_rel_diff_list:
    print("\n===== Sorted Isotopes by Average Relative Difference =====")
    for isotope, avg_rel in sorted(avg_rel_diff_list, key=lambda x: x[1]):
        print(f"{isotope:<10} {avg_rel:10.2f} %")