# Nuclear data



### Config the nuclear data

This minimal example uses a single nuclear data library for all nuclides.
Other configurations are covered under the nuclear data options below.

```python
yamc.cross_section_data = "endf-b8.1"
```

### Nuclear data settings

Five module-level settings, read and written like attributes. Each takes a library
keyword or a path, and each can point somewhere different, so a run can mix
libraries by subsection.

| setting | supplies | decides | left unset |
| --- | --- | --- | --- |
| ``cross_section_data`` | continuous-energy cross sections per nuclide | every reaction rate, in transport and in transmutation alike | nothing at all: there is no data to transport with, and every transmutation ``sigma * phi`` is zero |
| ``transmutation_reactions`` | which reaction on which nuclide gives which product, and its Q | which activation products exist to be made | ``endf-b8.1`` |
| ``transmutation_decay_data`` | half-lives, decay modes, mean decay energies, decay photon lines | how the inventory decays, and the photon source a shutdown dose rate transports | ``endf-b8.1`` |
| ``transmutation_fission_yields`` | fission product yields | what a fission makes | ``endf-b8.1`` |
| ``transmutation_branch_ratios`` | isomeric branching against incident energy | which product is left in a metastable state | no overlay, so every product lands in its ground state |

Four of the five have a default and give an answer without being set. ``cross_section_data`` has none, and every
transmutation and shutdown-dose example on these pages assumes you have set it.

Photons are chosen separately and default to ``endf-b8.1`` whatever the neutron
library is. See
[Nuclear data libraries](nuclear_data_libraries.md#photon-data-is-chosen-separately).

#### The two that both say reactions

``cross_section_data`` and ``transmutation_reactions`` are halves of one
calculation rather than alternatives. For each nuclide the network lists the
channels it has, each channel's name is mapped to an MT number, and the cross
sections are what that MT is priced with. So the network decides which channels
are considered at all and the cross sections decide how fast each one runs.
Swap the library behind ``transmutation_reactions`` and products appear or
disappear; swap the one behind ``cross_section_data`` and the same products
arrive at different rates.

They also fail differently, which is how to tell them apart in practice. A missing
``cross_section_data`` is silent: the inventory comes back as it went in. A
network pointed at the wrong material is refused before the solve, naming what
the material holds against what the network drives.

### Nuclide properties and cross sections

Nuclides can be created and their properties inspected.

```python
nuclide = yamc.Nuclide('Li6')

print(nuclide.name)                    # "Li6"
print(nuclide.element)                 # "Li"
print(nuclide.atomic_number)           # 3
print(nuclide.mass_number)             # 6
print(nuclide.available_temperatures)  # ["250", "294", "600", "900", "1200", "2500"]
print(nuclide.reaction_mts)            # [1, 2, 3, 4, 51, ..., 105]
```

Microscopic cross sections are retrieved by reaction name or MT number. Each call returns ``(xs_values, energy_grid)``. Different reactions may have different energy grids.

```python
xs_total, energy = nuclide.microscopic_cross_section(reaction=1, temperature="294")           # total
xs_elastic, energy = nuclide.microscopic_cross_section(reaction=2, temperature="294")         # elastic
xs_ngamma, energy = nuclide.microscopic_cross_section(reaction="(n,gamma)", temperature="294")
xs_nt, energy = nuclide.microscopic_cross_section(reaction="(n,t)", temperature="294")
```

Plotting cross sections with matplotlib (the figure below is rendered live from this code when the docs are built):

```python exec="1" html="1" source="above"
import matplotlib.pyplot as plt
import yamc

yamc.cross_section_data = "endf-b8.1"

nuclide = yamc.Nuclide("Li6")
reactions = [
    (1, "Total (MT 1)"),
    (2, "Elastic (MT 2)"),
    ("(n,gamma)", "(n,gamma)"),
    ("(n,t)", "(n,t)"),
]

fig, ax = plt.subplots(figsize=(8, 5))
for reaction, label in reactions:
    xs, energy = nuclide.microscopic_cross_section(reaction=reaction, temperature="294")
    ax.loglog(energy, xs, label=label)

ax.set_xlabel("Energy [eV]")
ax.set_ylabel("Microscopic Cross Section [barns]")
ax.set_title("Li6 microscopic cross sections")
ax.legend()
ax.grid(True, alpha=0.3)

from io import StringIO as _S
_b = _S(); plt.savefig(_b, format="svg", bbox_inches="tight"); print(_b.getvalue())
```

### Elements and natural abundances

Element data is available via the ``yamc.data`` module.

```python
# Natural abundances
abundances = yamc.data.natural_abundance()
print(abundances["Li6"])  # 0.07589
print(abundances["Li7"])  # 0.92411

# List nuclides for an element
nuclides = yamc.data.element_nuclides()
print(nuclides["Li"])  # ["Li6", "Li7"]
```

The abundances are the IUPAC TICE 2013 isotopic compositions (Meija et al.,
Pure Appl. Chem. 88 (2016) 293, doi:10.1515/pac-2015-0503). Each value is the
column 9 representative isotopic composition. The 12 elements that TICE gives
as an interval rather than a value (H, Li, B, C, N, O, Mg, Si, S, Cl, Br, Tl)
use the column 6 best measurement instead, and mononuclidic elements are 1.

``natural_abundance_records()`` returns the full TICE row behind each value:
the column 9 value and its expanded uncertainty (or the interval), the column 4
observed range of natural variation, and the column 6 best measurement with its
uncertainty, coverage factor and calibration flag. A field TICE leaves empty is
``None``, which means "not stated", never zero.

```python
records = yamc.data.natural_abundance_records()
fe58 = records["Fe58"]
print(fe58["representative_value"], fe58["representative_uncertainty"])  # 0.00282 0.00012
print(fe58["best_measurement"], fe58["best_measurement_uncertainty"])    # 0.002819 2.7e-05
print(fe58["best_measurement_coverage"])                                 # '2s'
print(records["Li6"]["representative_interval"])                         # (0.019, 0.078)
```

These uncertainties are reference data only. No calculation samples or
propagates them: turning them into a distribution needs choices TICE does not
make (a coverage factor for column 9, correlations between the isotopes of an
element, how to treat an interval).

The element symbol table is bound too, so a script naming files or directories
from ENDF data does not need its own copy of the periodic table.

```python
yamc.data.atomic_symbol(26)         # 'Fe'
yamc.data.atomic_number("Fe")       # 26
yamc.data.split_nuclide("Fe56")     # ('Fe', 56, 0)
yamc.data.split_nuclide("Ta180_m1") # ('Ta', 180, 1)
```

``split_nuclide`` returns ``(symbol, mass number, isomeric state)``, the state
being 0 for a ground state and 1 for the first metastable. Splitting a name by
stripping its trailing digits does not do this: that leaves ``"Ta180_m"`` with a
mass number of 1, and natural tantalum really does contain Ta180m, so a
composition built from ``natural_abundance()`` will hand you that name.

``atomic_symbol`` takes Z = 0 and answers ``'n'`` for the neutron rather than
raising, because ENDF libraries really do ship a neutron material (JENDL is one)
and a caller walking a sublibrary will meet it. ``atomic_number`` is case
sensitive and takes the symbol alone: ``"fe"`` and ``"iron"`` are both refused.

### Options for setting nuclear data source

The source of nuclear data can be set at the nuclide, material and program level in a few ways.

#### cross_section_data with a library keyword

The simplest method is to use a single source of nuclear data for everything.

```python
yamc.cross_section_data = "fendl-3.2d"
```

There are six recognized keywords: ``endf-b8.1``, ``jeff-4.0``, ``jendl-5.0``,
``tendl-2025``, ``tendl-2017`` and ``fendl-3.2d``.
[Nuclear data libraries](nuclear_data_libraries.md) has what each one publishes,
how big each one is, and which can serve as a photon source.

TENDL is neutron-only, so it has no decay or fission-yield evaluations; pair it with a photon library (``endf-b8.1`` by default) for coupled neutron/photon transport, and take ``decay``/``fission_yields`` from a library that has them while sourcing ``reactions`` and/or ``branching`` from a TENDL release. Pointing a transmutation part at a library that does not provide it (e.g. ``decay`` from a TENDL release) fails immediately with a message listing the subsections that library does provide.

Whenever nuclear data is needed, yamc checks the user's ```.cache/yamc``` folder for the Arrow data of the required nuclide.
If the file is there it is used; otherwise it is downloaded to the cache first.

#### cross_section_data with a directory path

If all your Arrow data directories live in one directory you can point `cross_section_data` at that directory.
When a nuclide is needed (e.g. `Li6`), YAMC resolves it to `<dir>/Li6.arrow` automatically.

<!-- doctest: skip -->
```python
yamc.cross_section_data = "/path/to/neutron/"
```

#### cross_section_data with file paths

Data files can also be downloaded from the nuclear data repository [FENDL-3.2d](https://github.com/fusion-neutronics/cross_section_data_fendl_3.2d_arrow). Once they are saved locally, their paths can be used to configure the nuclear data. Replace the `tests` dir with the path to the downloaded data files.

```python
yamc.cross_section_data = {
    "Be9": "tests/Be9.arrow",
    "Fe54": "tests/Fe54.arrow",
    "Fe56": "tests/Fe56.arrow",
    "Fe57": "tests/Fe57.arrow",
    "Fe58": "tests/Fe58.arrow",
    "Li6": "tests/Li6.arrow",
    "Li7": "tests/Li7.arrow",
}
```

#### cross_section_data with paths and specific libraries

Different sources can also be mixed. In this example we use the `endf-b8.1` and `fendl-3.2d` library keywords for some nuclides and local file paths for others. Replace the ```tests``` dir with the path to the downloaded data files.

```python
yamc.cross_section_data = {
    "Be9": "fendl-3.2d",
    "Fe54": "endf-b8.1",
    "Fe56": "tests/Fe56.arrow",
    "Fe57": "tests/Fe57.arrow",
    "Fe58": "tests/Fe58.arrow",
    "Li6": "endf-b8.1",
    "Li7": "fendl-3.2d",
}
```

#### Specific nuclear data on the Nuclide

You can also avoid setting ```yamc.cross_section_data``` entirely and specify the nuclear data to use on the Nuclide object directly.

This can be done using a nuclear data library.
```python
nuclide = yamc.Nuclide('Li6')
nuclide.read_nuclear_data('endf-b8.1')
```

Alternatively a specific nuclear data path

```python
nuclide = yamc.Nuclide('Li6')
nuclide.read_nuclear_data('tests/Li6.arrow')
```

Or a directory path (resolves to `tests/Li6.arrow`)

```python
nuclide = yamc.Nuclide('Li6')
nuclide.read_nuclear_data('tests')
```

#### Sample interacting reaction

Once the interacting nuclide is known then the reaction can be sampled using ```Nuclide.sample_reaction()```.

```python
import yamc
nuc = yamc.Nuclide('Li6')
nuc.read_nuclear_data('tests/Li6.arrow')
reaction = nuc.sample_reaction(energy=1.0, temperature='300', seed=42)
print(f'interacting reaction is {reaction}')
```

### Transmutation data

For transmutation (``Material.transmute()``), coupled transport-transmutation (``Model.simulate_transmutation()``) and shutdown dose rate (D1S), YAMC also needs a transmutation network: the chain of reactions, decays and decay-photon sources used to evolve nuclide inventories.

The network is assembled from four independently sourced parts, each set globally to a library keyword (its subsection tarball is downloaded and cached on first use) or a path to a local chain directory:

<!-- doctest: skip -->
```python
# decay / reactions / fission_yields default to the endf-b8.1 library when
# unset; set any part independently (keyword or path). Parts may come from
# different libraries.
yamc.transmutation_decay_data     = "endf-b8.1"
yamc.transmutation_reactions      = "endf-b8.1"
yamc.transmutation_fission_yields = "endf-b8.1"

# Isomeric branching is a separate overlay part; sourced (when available) from a
# TENDL release. Unset means no branching overlay is applied.
yamc.transmutation_branch_ratios  = "tendl-2025"

# The reactions network topology + Q values can also come from a TENDL release,
# which covers far more parent nuclides than endf-b8.1:
yamc.transmutation_reactions      = "tendl-2017"

# A value can also be a path to a local chain directory that holds decay/,
# reactions/ and fission_yields/ subsection subdirectories:
yamc.transmutation_decay_data = "/path/to/transmutation-endf-b8.1-sfr.arrow"
```

See [Nuclear data libraries](nuclear_data_libraries.md) for which library provides which subsection. Setting a part to a library that does not publish it fails fast with a message listing that library's available subsections.

A part can also be turned off, which is the correct setting for a decay-only
calculation or for a material in which nothing fissions:

<!-- doctest: skip -->
```python
yamc.set_transmutation_reactions(False)       # decay only, no reaction topology
yamc.set_transmutation_fission_yields(False)  # nothing here fissions
```

``False`` is not the same as leaving a part unset, which falls back to the default
library. A rate that then needs the part that was switched off is refused when the
burnup matrix is built, naming the nuclide and the reaction, rather than being
solved as though the reaction produced nothing, or losing its fission products
while still burning the nuclide.

The transmutation APIs and the decay-photon post-processing (``PulseSchedule.time_correct_tally``) read these settings; there is no per-call override.

