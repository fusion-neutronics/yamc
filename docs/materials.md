# Materials

### Creating a material

A material is defined by a composition dictionary and a density. The simplest form uses nuclide names as keys with atom fractions as values.

```python
iron = yamc.Material(
    composition={
        "Fe54": 0.0585,
        "Fe56": 0.9175,
        "Fe57": 0.0212,
        "Fe58": 0.0028,
    },
    density=7.874,
)
```

### Composition keys

Composition keys can be nuclide names, element symbols, or chemical formulas. Elements are automatically expanded to their natural isotopic distribution and formulas are parsed into their constituent nuclides.

```python
# Nuclide names (explicit isotopes)
mat = yamc.Material(
    composition={"Li6": 0.5, "Li7": 0.5},
    density=0.534,
)

# Element symbols (expanded to natural isotopes)
mat = yamc.Material(
    composition={"Fe": 1.0},
    density=7.874,
)

# Chemical formulas
mat = yamc.Material(
    composition={"H2O": 1.0},
    density=1.0,
)

# Mixed keys
mat = yamc.Material(
    composition={"Li4SiO4": 0.6, "Be": 0.4},
    density=2.0,
)
```

### Density specification

Density is always a float; its meaning is set by the optional ``units`` keyword. With ``units="atom/barn-cm"`` the value is the **total atom density** and the composition values are relative atom fractions (requires ``percent_type="atom"``).

```python
# Mass density (default units g/cm3)
mat = yamc.Material(composition={"Fe": 1.0}, density=7.874)
mat = yamc.Material(composition={"Fe": 1.0}, density=7874.0, units="kg/m3")

# Total atom density + relative atom fractions
mat = yamc.Material(composition={"Li6": 0.625, "Li7": 0.375}, density=0.08,
                    units="atom/barn-cm")

# Absolute per-nuclide atom densities, stored verbatim (e.g. from a solver)
mat = yamc.Material.from_atom_densities({"Li6": 0.05, "Li7": 0.03})
```

``material.density`` is always a float -- the mass density for ``g/cm3``/``kg/m3``, or the total atom density for ``atom/barn-cm`` (and for transmuted materials).

### Atom vs weight percent

The ``fraction_type`` parameter controls how composition fractions are interpreted. The default is ``"atom"`` (atom fractions).

```python
# Atom fractions (default)
mat = yamc.Material(
    composition={"U235": 0.05, "U238": 0.95},
    density=19.1,
    fraction_type="atom",
)

# Mass fractions
mat = yamc.Material(
    composition={"U235": 0.05, "U238": 0.95},
    density=19.1,
    fraction_type="mass",
)
```

### Enrichment

The ``enriched()`` helper lets you specify isotopic enrichment for elements and formulas.

```python
# Lithium enriched to 60% Li6
mat = yamc.Material(
    composition={
        "Li": yamc.enriched(1.0, target="Li6", percent=60.0),
    },
    density=0.534,
)

# Mixed: enriched element + natural element
mat = yamc.Material(
    composition={
        "Li": yamc.enriched(0.4, target="Li6", percent=90.0),
        "Fe": 0.6,
    },
    density=5.0,
)

# Mass percent enrichment
mat = yamc.Material(
    composition={
        "Li": yamc.enriched(1.0, target="Li6", percent=90.0, fraction_type="mass"),
    },
    density=0.534,
)
```

### Material properties

The constructor accepts several optional keyword arguments for naming, identification, temperature selection, and transmutation.

```python
mat = yamc.Material(
    composition={"Fe": 1.0},
    density=7.874,
    name="iron",
    temperature=294,    # Kelvin - selects temperature in nuclear data
    transmutable=True,   # required for Model.simulate_transmutation()
    volume=523.6,       # cm^3 - required for transmutation normalization
)
```

Properties are readable after construction.

```python
print(mat.name)                  # "iron"
print(mat.nuclides)             # [("Fe54", 0.0585), ("Fe56", 0.9175), ...]
print(mat.density)              # 7.874
print(mat.density_units)        # "g/cm3"
print(mat.temperature)          # "294"
print(mat.average_molar_mass()) # g/mol
```

### Transmutable materials

Only materials marked ``transmutable=True`` are evolved by ``Model.simulate_transmutation()``. Inert materials (the default) keep their composition fixed across timesteps, which avoids needlessly tracking nuclides in shielding or structural cells. A ``volume`` (cm³) is also required so reaction rates can be normalised to atom densities.

```python
fuel = yamc.Material(
    composition={"U235": 0.05, "U238": 0.95},
    density=10.4,
    name="fuel",
    transmutable=True,   # evolved during transmutation
    volume=523.6,        # cm^3 - required when transmutable=True
)

shielding = yamc.Material(
    composition={"Fe": 1.0},
    density=7.874,
    name="shielding",
    # transmutable=False is the default -- composition stays fixed
)
```

The flag is not needed for ``Material.transmute()`` (standalone transmutation with a pre-computed flux) -- that method always evolves the material it is called on.

### Mixing materials

Existing materials can be combined with ``Material.mix_materials()``. The ``fraction`` parameter controls whether fractions are by atom (``"atom"``), mass (``"mass"``), or volume (``"volume"``).

```python
steel = yamc.Material(composition={"Fe": 1.0}, density=7.874)
water = yamc.Material(composition={"H2O": 1.0}, density=1.0)

# 80% steel, 20% water by volume
mixed = yamc.Material.mix_materials(
    materials=[steel, water],
    fractions=[0.8, 0.2],
    fraction_type="volume",
    name="steel_water_mix",
)
```

### Material collections

Rather than looking up a composition and typing it in, a material can be pulled straight from a curated collection. ``yamc.materials.pnnl`` holds 410 materials from the [PNNL Compendium of Material Composition Data for Radiation Transport Modeling (PNNL-15870, Rev. 2)](https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev2.pdf): concretes, steels, tissues, air, soils, shielding, scintillators, fuels and more. The data is bundled in the wheel, so a lookup needs no network and cannot fail.

Materials are keyed by the report's own name, verbatim, so an entry can be quoted straight back to it.

```python
steel = yamc.materials.pnnl["Steel, Stainless 304"]
```

The result is an ordinary ``Material``, built fresh on every lookup, so it can be renamed and mutated freely. Renaming is what you want when a mesh geometry expects a particular material name:

```python
steel = yamc.materials.pnnl["Steel, Stainless 304"]
steel.name = "firstwall_material"
steel.temperature = 600
```

See the [full list of PNNL materials](pnnl_materials.md) for every name, or search the collection from Python. ``search()`` is a case-insensitive substring match and returns names you can pass straight back in.

```python
yamc.materials.pnnl.search("concrete")   # -> ['Concrete [Los Alamos (MCNP) Mix]', ...]
yamc.materials.pnnl.names()              # -> all 410 names
len(yamc.materials.pnnl)                 # -> 410
"Steel, Stainless 304" in yamc.materials.pnnl
yamc.materials.collections()             # -> ['pnnl']
```

An unrecognised name raises ``KeyError`` suggesting near matches rather than failing silently.

```python
yamc.materials.pnnl["Steel, Stainless 999"]
# KeyError: no material named 'Steel, Stainless 999' in the pnnl collection;
#           did you mean: Steel, Boron Stainless, Steel, HT9 Stainless, ...
```

#### Overriding fields at lookup

``temperature``, ``volume``, ``transmutable``, ``name`` and ``id`` are all settable on the returned material. ``density`` is not, because a ``Material``'s per-nuclide atom densities are fixed at construction, so use ``material()`` to build at a different density, e.g. a concrete at a different porosity.

```python
shield = yamc.materials.pnnl.material(
    "Concrete, Ordinary (NIST)",
    name="bioshield",
    density=1.8,        # g/cm3, replacing the collection's reference density
    temperature=600,
    volume=1.0e6,
    transmutable=True,
)
```

#### Provenance

Every entry carries the report's own entry number so a composition can be cited.

```python
yamc.materials.pnnl.entry("Steel, Stainless 304")
# {'name': 'Steel, Stainless 304', 'number': 331, 'density': 8.03,
#  'formula': None, 'composition': {'C': 0.003635, 'Mn': 0.01987, ...},
#  'citation': 'PNNL-15870, Rev. 2',
#  'url': 'https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev2.pdf'}
```

#### How compositions are stored

Values are the report's own **atom fractions**, taken verbatim from its Elemental column, so no weight-to-atom conversion is applied on the way in.

Where an element is isotopically natural the **element symbol** is stored, and yamc expands it with yamc's natural abundances rather than pinning the report's. Where a material is genuinely not natural (enriched lithium and boron, uranium, plutonium, heavy water, He-3) the **individual nuclides** are stored instead.

```python
yamc.materials.pnnl.entry("Steel, Stainless 304")["composition"]  # {'Fe': 0.667971, 'Cr': 0.199443, ...}
yamc.materials.pnnl.entry("Water, Heavy")["composition"]          # {'H2': 0.666667, 'O': 0.333333}
```

!!! note "One entry is omitted"
    The report has 411 entries; yamc bundles 410. Entry 173 "Iron Boride (Fe2B)" is omitted because it contradicts itself: its stated formula and molecular weight are those of Fe₂B, while its atom fractions, weight fractions and atom densities all describe FeB₂. Entry 174 "Iron Boride (FeB)" is unaffected and is included.

The compositions are reproduced from a US Department of Energy technical report (R.J. McConn Jr et al., Pacific Northwest National Laboratory, 2021). yamc bundles them for lookup convenience and claims no authorship of them.

### Loading nuclear data

Nuclear data can be loaded automatically via ``yamc.cross_section_data`` or explicitly per material.

<!-- doctest: skip -->
```python
# Automatic: yamc.cross_section_data downloads and caches data as needed
yamc.cross_section_data = "fendl-3.2d"

# Explicit: per-material file paths
mat.read_nuclear_data({
    "Fe54": "path/to/Fe54.arrow",
    "Fe56": "path/to/Fe56.arrow",
})

# With photon data (required when transport_secondary_photons=True on Model)
mat.read_nuclear_data(
    {"Fe56": "path/to/Fe56.arrow"},
    photon_data={"Fe": "path/to/Fe.arrow"},
)
```

Arrow (``.arrow``) format is used for both neutron and photon data. Photon data uses the element symbol as the key (not the nuclide name).

### Macroscopic cross sections

The macroscopic cross section for a reaction can be retrieved by name or MT number. Returns ``(xs_values, energy_grid)`` on the material's unified energy grid.

<!-- doctest: skip -->
```python
# By reaction name
xs, energy = mat.macroscopic_cross_section(reaction="(n,total)")

# By MT number
xs, energy = mat.macroscopic_cross_section(reaction=1)   # total
xs, energy = mat.macroscopic_cross_section(reaction=102) # (n,gamma)
```

Plotting macroscopic cross sections (rendered live when the docs are built):

```python exec="1" html="1" source="above"
import matplotlib.pyplot as plt
import yamc

yamc.cross_section_data = "endf-b8.1"

mat = yamc.Material(composition={"Fe": 1.0}, density=7.874, temperature=294, name="iron")

reactions = [(1, "Total"), (2, "Elastic"), ("(n,gamma)", "(n,gamma)")]
fig, ax = plt.subplots(figsize=(8, 5))
for reaction, label in reactions:
    xs, energy = mat.macroscopic_cross_section(reaction=reaction)
    ax.loglog(energy, xs, label=label)

ax.set_xlabel("Energy [eV]")
ax.set_ylabel("Macroscopic Cross Section [1/cm]")
ax.set_title("Iron macroscopic cross sections")
ax.legend()
ax.grid(True, alpha=0.3)

from io import StringIO as _S
_b = _S(); plt.savefig(_b, format="svg", bbox_inches="tight"); print(_b.getvalue())
```

### Mean free path

The neutron mean free path at a given energy (in eV) is returned in cm.

<!-- doctest: skip -->
```python
mfp = mat.mean_free_path_neutron(14e6)  # 14 MeV
```

### Sampling

Monte Carlo transport works by sampling random distances and interactions. These methods expose that process for inspection or educational use.

<!-- doctest: skip -->
```python
# Sample a distance to the next collision (cm)
distance = mat.sample_distance_to_collision(energy=14e6, seed=42)

# Sample which nuclide interacts
nuclide_name = mat.sample_interacting_nuclide(energy=14e6, seed=42)
```

### Transmutation

Materials can be transmuted standalone using a pre-computed multigroup flux spectrum, without running transport. This is useful when the flux is known from a prior simulation.

The example below uses a simple hand-specified 3-group spectrum so it can be read in isolation. In a real workflow the spectrum is usually obtained from a transport tally -- see [Neutron and photon flux spectra](tallies.md#neutron-and-photon-flux-spectra) for how to score an energy-binned flux on a cell and feed ``results[tally].mean`` (scaled by the source strength) into ``transmute()``.

<!-- doctest: skip -->
```python
# The transmutation network is assembled from the per-subsection sources;
# point the three network-forming parts at the chain file.
yamc.transmutation_decay_data = "transmutation-endf-b8.1-sfr.arrow"
yamc.transmutation_reactions = "transmutation-endf-b8.1-sfr.arrow"
yamc.transmutation_fission_yields = "transmutation-endf-b8.1-sfr.arrow"

# The spectrum rides on each irradiation pulse: a NeutronSource whose energy is
# a Histogram over the energy groups. The Histogram normalizes the shape, so a
# multigroup flux from a tally can be passed straight in -- see "Neutron and
# photon flux spectra" for the full workflow. The pulse rate is the total flux
# magnitude [n/cm^2/s]; Cooldown is decay-only.
energy_groups = [1e-5, 0.625, 1e5, 2e7]   # 3 groups (eV, ascending); or a name like "CCFE-709"
multigroup_flux = [1e12, 5e12, 1e14]      # thermal, epithermal, fast

spectrum = yamc.NeutronSource(energy=yamc.sources.Histogram(energy_groups, multigroup_flux))
total_flux = sum(multigroup_flux)
schedule = yamc.PulseSchedule(
    [yamc.Pulse(rate=total_flux, duration=(1, "h"), source=spectrum) for _ in range(5)]  # 5h irradiation
    + [yamc.Cooldown(duration=(1, "h"))]                                                  # 1h cooling
)
results = material.transmute(schedule=schedule)

# results is a TransmutationResults. step_materials() gives one Material per
# schedule step, in order.
for i, mat in enumerate(results.step_materials(material_id=material.id or 0)):
    print(f"Step {i+1}: {len(mat.nuclides)} nuclides")

# The end state on its own, without walking the steps.
final = results.get_final_material(material_id=material.id or 0)
if final is not None:
    print(f"Final: {len(final.nuclides)} nuclides")
```

Parameters:
- ``schedule``: a ``PulseSchedule`` of ``Pulse`` / ``Cooldown`` steps. Each irradiation ``Pulse`` carries a duration (seconds, or a ``(value, unit)`` tuple), a ``rate`` (the total flux magnitude [n/cm²/s]), and a ``source`` -- a ``NeutronSource`` whose energy is a ``Histogram`` giving the multigroup spectrum (its boundaries may be a named group structure such as ``"CCFE-709"``). The source's position/direction are ignored (no transport runs). ``Cooldown`` is decay-only. Different pulses may carry different spectra.
The transmutation network is assembled from the configured per-subsection sources (``yamc.transmutation_decay_data``, ``yamc.transmutation_reactions``, ``yamc.transmutation_fission_yields``, and optionally ``yamc.transmutation_branch_ratios``).

#### Repeated calls, and letting the data go

``transmute()`` loads the cross sections it needs **into the material** and keeps
them there, so calling it again on the same material does no file reading at all.
On a steel against ENDF/B-VIII.1 that is 1.8 s of a 2.3 s call. The composition is
not modified.

That is the right trade for a material transmuted more than once, and the wrong
one for a sweep over thousands of *distinct* compositions, where every material
would hold a few hundred nuclides' data for as long as it lives.
``release_nuclear_data()`` drops it. The material stays usable, and the next call
that needs the data loads it again.

<!-- doctest: skip -->
```python
for composition in many_compositions:
    scan_mat = yamc.Material(composition=composition, density=7.9)
    scan_mat.transmute(schedule=schedule)
    scan_mat.release_nuclear_data()
```

#### Resonance self-shielding

The spectrum handed to ``transmute()`` is the flux *outside* the material. Inside
a lump the flux is depressed at the energies where the material's own total cross
section is large, so a resonance absorber shields itself and its reaction rates
come out lower than a dilute fold predicts.

Give ``transmute()`` a lump size and that depression is computed and applied:

<!-- doctest: skip -->
```python
results = material.transmute(
    schedule=schedule,
    self_shielding_shape=yamc.shapes.FoilLump(thickness=0.1),   # cm
)
```

The five shapes live in ``yamc.shapes``. The shape decides how much surface a
volume hides behind, since the mean chord is ``4V/S``.
``SphereLump()`` and ``CubeLump()`` take no arguments because
``Material(volume=...)`` already fixes them; ``FoilLump(thickness=)``,
``CylinderLump(radius=)`` and ``WireLump(radius=)`` carry the one dimension a
volume cannot imply. A sphere has the least surface for its volume, so it has the
longest chord and shields more than any other shape of the same size: an upper
bound rather than a safe default, which is why there is no default at all.

``self_shielding_chord`` takes the mean chord ``4V/S`` in cm directly, for a lump
none of the five shapes describes. Give one or the other, not both.

Omit both and nothing is shielded: the rates are bit-identical to a run without
the argument, and no shape is inferred from a geometry this material does not
have.

The flux inside the lump comes from a slowing-down solve, which assumes nothing
about resonances being narrow. The cheaper narrow-resonance approximation is
deliberately not offered, because it over-shields strong elastic scatterers badly
enough to be worse than applying no correction at all. W186(n,gamma) in FNG-tung
measures 1.29 b, where the slowing-down solve gives C/E 0.80, no correction gives
0.86 and narrow resonance gives 0.46.

``get_self_shielding_info`` reports what was done, and is present for every
material a ``transmute()`` call solved, shielded or not:

<!-- doctest: skip -->
```python
shielding = results.get_self_shielding_info(material_id=material.id or 0)
if shielding is not None:
    shielding["method"]            # how the depression was computed
    shielding["chord_cm"]          # None if the run was dilute
    shielding["shielded"]          # nuclides the correction reached
    shielding["not_shielded"]      # the ones it did not, and why
    shielding["strongest_factor"]  # smallest factor any group average was scaled by
```

A ``strongest_factor`` of ``1.0`` means the run shielded nothing in practice,
which is a different statement from not having tried.

A **dilute** run fills ``would_shield`` instead: the nuclides whose own resonances
are structured enough to have suppressed a reaction, each mapped to how strongly.
It is an indicator, not a correction and not a bound. It is computed from the one
reaction and that nuclide's own density, with no geometry in it, so it can sit
either side of the real factor: a lump thinner than a centimetre of chord shields
less, and a material whose other nuclides dominate the total at the resonance dips
the flux further. On the FNS tungsten foil it reads 0.634 for W186, where the
slowing-down correction on the same foil and spectrum gives 0.730. So it says
"this answer may be high, and this is a resonance absorber", which is the warning
a dilute run should carry rather than silence. Read the size of the effect off a
shielded run.

#### Nuclear-data uncertainty

Pass a ``DataUncertainty`` and every nuclide density comes back with a standard
deviation beside it:

<!-- doctest: skip -->
```python
results = material.transmute(
    schedule=schedule,
    data_uncertainty=yamc.DataUncertainty(seed=42),
)

mid = material.id or 0
# Both answer None for a nuclide, material or step the run does not have, so a
# real script checks rather than formatting None into the output.
mean = results.get_nuclide_density(material_id=mid, nuclide="Mn56", step=1)
sigma = results.get_nuclide_uncertainty(material_id=mid, nuclide="Mn56", step=1)
if mean is not None and sigma is not None:
    print(f"Mn56: {mean:.3e} +/- {sigma:.3e}  ({100 * sigma / mean:.1f}%)")
```

Each sample redraws the inputs listed under [Choosing what to
perturb](#choosing-what-to-perturb) from their stated uncertainties: the
activation cross sections from their ENDF MF=33 covariance, folded against this
material's own spectrum; the half-lives and decay energies from the decay data's
own standard deviations; and the flux from the error you give it. The schedule is
re-solved with the drawn rates, flux and half-lives, and the drawn decay energies
are applied where decay heat is evaluated. That is exact to all orders in the
matrix exponential: nothing is linearized, and the solver itself is untouched,
only its input changes.

``DataUncertainty(seed=...)`` is reproducible. A nuclide's cross-section,
half-life and decay-energy draws in a given sample are a pure function of
``(seed, sample, nuclide)``, and the flux draw of ``(seed, sample, spectrum)``,
the spectrum's place among the schedule's spectra. The statistical draw is keyed
on ``(seed, sample)`` and made jointly over the material's tallied rates, so
adding or removing a tallied rate changes the draws of the others. In every case
the same seed gives the same answer whatever the sample count or iteration
order. Leave ``samples`` unset and the solver adds samples until the reported
sigmas settle; setting it bounds cost or reproduces a specific run, but it is
not an accuracy dial.

Omitting the argument costs nothing. No covariance is read from disk, no matrix
is folded or factorized, and the inventories are bit-identical.

##### Choosing what to perturb

Five sources can be perturbed, and by default every one that applies to the call
is:

<!-- doctest: skip -->
```python
yamc.DataUncertainty.available_sources()
# ['cross_sections', 'flux_spectrum', 'half_life', 'statistical', 'decay_energy']

# isolate one contribution
xs_only = yamc.DataUncertainty(seed=42, sources=["cross_sections"])
```

- ``cross_sections``: the activation cross sections, from the evaluation's MF=33
  covariance folded against the spectrum the rates are collapsed with.
- ``flux_spectrum``: the error on the spectrum you supplied, described below.
  The spectrum calls (``Material.transmute`` and ``yamc.transmute``) only. It
  moves every rate collapsed from the spectrum's groups, which leaves out the
  ``(n,n')`` channels the isomeric overlay adds (see [What it does not
  cover](#what-it-does-not-cover)).
- ``half_life``: each reachable unstable nuclide's half-life, from the decay
  data's standard deviation on it. A sample's half-lives are used in its solve
  and again in the activity, decay heat and dose evaluated from its inventory.
  That keeps a saturated activity (``lambda N = R``) insensitive to its own
  half-life, where using the nominal value in either place would hand it the
  density's whole spread. It is the only source a shutdown-dose
  ``PulseSchedule.time_correct_tally`` applies.
- ``statistical``: the Monte Carlo uncertainty of transport-tallied reaction
  rates, from their per-history covariance, so rates scored by the same
  histories move together.
  [Model.simulate_transmutation()](simulate_transmutation.md) only, which
  accepts ``data_uncertainty`` with ``method="independent"`` alone:
  ``method="coupled"`` raises.
- ``decay_energy``: each nuclide's mean decay energy, from the sigmas the decay
  data gives its recoverable-heat components (beta, gamma, alpha), or from the
  total's where the converted data carries only the total. The evaluation
  states no correlation between the components, so they are drawn as
  independent and the total's sigma is taken as their quadrature sum. A decay
  energy never enters the solve, so this moves decay heat and nothing else.

A call ignores the sources that do not apply to it, and ``sources`` in the report
below lists the ones that did.

Restricting ``sources`` is how a run isolates one contribution, so that adding a
source and watching the inventory sigma grow is a measurement rather than a guess.
Naming a source this build cannot perturb raises rather than being ignored: a
source that has not landed yet must not look like one that contributed nothing.

``flux_spectrum`` propagates the error bar on the spectrum you supplied, which
means the pulse has to carry one. ``Pulse(flux_std_dev=...)`` takes the per-bin
standard deviation, in the same units and the same order as the ``Histogram``
values; from a yamc run that is ``results[tally].standard_deviation``, scaled by
whatever factor the mean was. The error is taken relative to the ``Histogram``
values as given, so a mean multiplied by the source strength beside an unscaled
standard deviation reads a flux error smaller by that factor, and nothing flags
it.

<!-- doctest: skip -->
```python
pulse = yamc.Pulse(
    rate=total_flux,
    duration=(1, "h"),
    source=spectrum,
    flux_std_dev=flux_sigma,   # per bin, same order as the Histogram values
)
```

A per-bin sigma treats the bins as independent, and the bins of a Monte Carlo
tally are not: the same histories score them, so they move together, and their
standard deviations alone understate the error of every rate that sums over a
band. ``Pulse(flux_covariance=...)`` takes the full covariance of the
``Histogram`` values instead, in the square of their units, which a tally built
with ``Tally(covariance=True)`` reports as ``results[tally].covariance``. It
scales with the square of the factor the mean was scaled by: a mean multiplied by
the source strength ``S`` needs the covariance multiplied by ``S**2``. Give one or
the other.

Omitting both is the common case: a spectrum taken from a published reference set
carries no stated error. The ``flux_spectrum`` source then contributes nothing and
says so in the report's ``spectra_without_flux_sigma``, rather than reading as a
flux known exactly.

##### How a sample is drawn

A cross-section channel's rate is multiplied by a lognormal factor with mean one
and that channel's variance in the matrix it samples from (the folded covariance,
or the repaired one described below), driven by the channel's correlated Gaussian
deviate. So each channel keeps that mean and variance exactly and cannot go
negative. Between channels the draw keeps the rank correlation; the
linear correlation shifts, by an amount that vanishes as the sigmas go to zero
and grows with them.

The other sources draw from a Gaussian on the stated sigma or covariance. A flux
bin, tallied rate or decay energy drawn below zero is set to zero, and a
half-life drawn at or below zero is set to a millionth of its nominal value.
``get_data_uncertainty_info`` counts the flux, tallied-rate and half-life cases
(``flux_bins_floored``, ``statistical_floored`` and ``half_lives_floored``, each
beside its ``_sampled`` total). A count well above zero says the sigma is a large
enough share of the value that a Gaussian no longer describes it. A shutdown-dose
``time_correct_tally`` floors its half-life draws the same way and does not count
them.

Evaluated covariance matrices are not always positive semi-definite, and one that
is not has no square root to sample from. The negative eigenvalues of the folded
matrix are set to zero, which is the standard repair, and the repair adds
variance rather than removing it: no channel's sampled sigma falls below the
evaluated one, and a channel the evaluation gave no variance can pick some up.
``matrices_clipped`` counts the matrices repaired and ``worst_relative_clip``
gives the largest ``|lambda_min| / lambda_max`` among them. That ratio is set by
the widest channel, often one carrying almost none of the rate, so a small value
does not mean the channel that matters kept its evaluated sigma.

##### What it does not cover

Every input below is held at its evaluated or nominal value in every sample, or
drawn without a correlation the data states. Some carry published uncertainties
the converted data does not hold yet, some have none to read, and some are a
property of the model rather than of the data. Most leave a term out of the
reported spread, so the sigma is too small by it. Two can go either way, because
holding them also removes a feedback that narrows the spread: self-shielding,
and the flux's response to the data. Their items say which way each goes, and
where.

- **Decay branching ratios.** The decay data states a sigma on many of them, and
  it is not in the converted data yet.
- **Fission yields.** The evaluations state a sigma on each yield, and it is not
  in the converted data yet.
- **Isomeric branching.** The split between a product's ground and metastable
  states is held. The evaluations that publish radionuclide-production
  covariance (MF=40) give the covariance of each state's production cross
  section, which is parsed but not yet in the converted data ([Covariance is
  presence, not
  coverage](nuclear_data_libraries.md#covariance-is-presence-not-coverage)).
  ENDF/B-VIII.1, FENDL-3.2d and TENDL state no correlation there between a
  product's states, or with the reaction's MF=33, so even once it is stored the
  split's uncertainty is only partly stated. With
  ``transmutation_branch_ratios`` configured, the ``(n,n')`` channels the
  overlay adds to reach a metastable state, In115 to In115m among them, come
  from its MF=10 partials with no group-averaged cross section behind them. They
  draw no cross-section covariance, and on the spectrum calls no
  ``flux_spectrum`` perturbation either, so they are missing from
  ``rate_fraction_covered``.
- **Resonance-parameter covariance (MF=32).** It is not read. Many evaluations
  keep the resolved-resonance uncertainty on capture, elastic and fission there
  and leave MF=33 at zero over that range, which reaches down to thermal. A
  product made by capture at thermal or resonance energies then draws little or
  none of its cross section's uncertainty. ENDF/B-VIII.1 W186 is one: its
  resonance parameters carry about 1.5% on capture below 10 keV, where its MF=33
  gives zero. Where the MF=33 block is stated as zero over that range, as
  W186's is, the report does not flag it yet: a block stated as zero counts as
  covered, so the nuclide is listed under ``perturbed`` and that channel's
  ``rate_fraction_covered`` reads 1.
- **MF=33 blocks that are not used.** Four kinds are left out. Blocks that
  correlate with another evaluation (``MAT1 != 0``, the links to the standards
  among them) are skipped, so each nuclide's cross sections are drawn
  independently of every other nuclide's; the same test also skips some blocks
  whose ``MAT1`` names the evaluation itself, which loses correlations within one
  nuclide. Blocks stated through other reactions' covariances (NC) are skipped,
  so a reaction whose covariance is given only that way draws none. Lumped
  reactions (MT=851 to 870), which state the uncertainty of a sum of reactions,
  are skipped too, so their components draw none of it. Tungsten's ``(n,2n)``,
  the channel behind most of a tungsten first wall's decay heat, is one:
  ENDF/B-VIII.1, JEFF-4.0 and FENDL-3.2d state its covariance only lumped with
  ``(n,2np)``. Blocks stated on partial levels (MT=600 to 849 and 875 to 891)
  are not used either, since the chain drives the total, so a reaction whose
  covariance is given only per level draws none: ENDF/B-VIII.1 Ca40 states its
  ``(n,p)`` and ``(n,a)`` only on MT=600 and 800, and is listed under
  ``perturbed`` through its capture alone. ``skipped_cross_material`` and
  ``skipped_nc`` count the first two kinds, per spectrum folded and including
  blocks for reactions the run does not drive. The lumps and the partial levels
  are not counted, and show only as a channel missing from
  ``rate_fraction_covered`` and a lower ``rate_fraction_covered_total``.
- **MF=33 blocks read in part.** The ENDF parser splits an LB=0 to 4 block's
  table in the wrong place, so an LB=0 to 2 block is used over only the lower
  half of its energy grid, and a channel whose uncertainty is stated above that
  draws little or none of it. 773 ENDF/B-VIII.1 blocks are affected: at 14 MeV
  Ni58 ``(n,a)`` draws 0.0% where the evaluation gives 19.8%, and Cr52 ``(n,p)``
  0.4% against 17.1%. The LB=4 expansion also swaps its two tables, which
  affects one block, FENDL-3.2d Ni58 ``(n,p)``, whose 14 MeV sigma reads 20.6%
  against the evaluation's 18.3%. The report does not flag either:
  ``rate_fraction_covered`` is taken from the grid as read, so the dropped part
  shows at most as rate outside the evaluation's grid, and not at all where
  another block on the same reaction spans the whole range, as one of Cr52's
  ``(n,p)`` blocks does.
- **Self-shielding.** With ``self_shielding_shape`` or ``self_shielding_chord``
  set, every sample uses the flux depression computed from the evaluated cross
  sections, and that pulls two ways. A perturbed capture cross section neither
  deepens nor relieves its own dip, and that feedback narrows the spread, since
  a larger capture then absorbs from a smaller flux. The elastic and total
  covariance never reaches the depression either, and that term widens the
  spread where the evaluation states a resonance-range elastic uncertainty.
  Which wins depends on the evaluation. For Au197 capture in a 0.1 mm foil under
  a 1/E spectrum on VITAMIN-J-175, ENDF/B-VIII.1 states no such elastic
  uncertainty, and the held sigma is too large: 1.70% held against 0.96% with
  the feedback. TENDL-2017 states one, and there the held sigma is too small:
  5.34% held against 6.07% with the feedback and the elastic term together.
  Those figures fold the shielded partial rates. The covariance fold does not
  yet: it uses the dilute partial rates while the rate itself is shielded, which
  overstates the relative sigma of a channel that shields strongly, on top of
  the above: 4.61% against the 1.70% for the ENDF/B-VIII.1 foil.
- **The flux's response to the data.** A perturbed cross section does not change
  the flux it is folded against. In ``Material.transmute`` the flux is your
  input, and in ``Model.simulate_transmutation`` the covariance is folded against
  the one tallied spectrum without re-running transport per sample. The term
  left out is the flux's first-order response to each perturbed channel, which
  grows with how much that channel shapes the flux. It is negligible for a trace
  nuclide or a weak channel, and not for a material that shapes its own flux or
  for the main constituents of a shield. A larger Li6(n,t) cross section in a
  breeder depresses the flux it sees, the same narrowing feedback as
  self-shielding, so holding the flux overstates that rate's sigma. A shield's
  cross sections move the flux of every product behind it, however trace, and
  leaving that out understates those products' sigmas. Elastic scattering, which
  shapes the flux without transmuting anything, is not drawn at all.
- **Decay photon line intensities.** Each line's emission probability per decay
  is held; a sample's half-life still scales its emission rate. The decay data
  states a sigma on nearly every line, and it is not in the converted data yet,
  so the nuclear-data spread on the decay-photon spectrum, contact dose and a
  shutdown dose rate carries nothing from the lines themselves.
- **Dose constants.** Contact dose is built from the material's photon mass
  attenuation coefficients, the air mass energy-absorption coefficients or the
  ICRP-116 conversion coefficients (by ``dose_quantity``), and a constant
  build-up factor (``build_up``, 2.0 by default). A shutdown dose tally is scored
  with fluence-to-dose coefficients. None of these tables publishes a per-value
  uncertainty and the build-up factor is a model constant, so none is perturbed.
- **Natural isotopic abundances.** An element is expanded into its isotopes with
  the same abundances in every sample. The published isotopic-composition tables
  state an uncertainty or a range of natural variation on every isotope of a
  multi-isotope element, above 1% of the value for about a third of them (Li6's
  stated range spans a factor of four). The bundled table carries one value per
  isotope and no uncertainty or range: the tables' representative value, or for
  the twelve elements they give as an interval, a single best measurement, which
  for Li6 (0.07589) sits near the top of its range.
- **Element fractions and density.** These are your input and start every sample
  at the values you gave. There is no argument to give their uncertainty, so it
  is not perturbed, and propagating it means sampling compositions yourself, one
  run per composition. For dose from a trace impurity it can outweigh everything
  above: in a 316L-type steel with 0.1 wt% cobalt, irradiated for five years at
  1e14 n/cm^2/s in an illustrative first-wall spectrum, the contact dose at ten
  years' cooling is almost all Co60, about 70% of it made from the cobalt itself,
  so a 30% uncertainty on the cobalt content is about 21% on the dose.

The report's ``not_perturbed`` names decay branching, fission yields, isomeric
branching and cross-material covariance (``MAT1 != 0``) in every run, and
half-life and decay energy when ``sources`` leaves them out. The rest of this list
applies whether or not the report names it.

Because a reported sigma of zero can mean "well known" or "nothing published",
the report, ``get_data_uncertainty_info``, keeps a nuclide with no usable MF=33
(``no_covariance_data``: none in the evaluation, or none this run could use,
such as blocks all skipped as above or all for reactions the run does not drive)
apart from one that states it. It does not yet separate a channel whose stated
covariance is zero, as in the MF=32 case above:

<!-- doctest: skip -->
```python
info = results.get_data_uncertainty_info(material_id=mid)
if info is not None:            # None unless data_uncertainty was passed
    info["sources"]             # the sources that applied to this run
    info["perturbed"]           # nuclides that had usable MF=33 covariance
    info["no_covariance_data"]  # nuclides with none usable, absent or skipped
    info["rate_fraction_covered"]  # per channel, the share of the rate the
                                   # covariance grid spans, zero blocks included
    info["rate_fraction_covered_total"]  # ... over the run, by production
    info["skipped_cross_material"]    # MF=33 blocks with MAT1 != 0, not used
    info["skipped_nc"]                # MF=33 blocks stated as NC, not used
    info["half_lives_perturbed"]      # nuclides whose half-life was drawn
    info["no_half_life_uncertainty"]  # ... and those whose data states no sigma
    info["decay_energies_perturbed"]  # nuclides whose decay energy was drawn
    info["no_decay_energy_uncertainty"]  # ... and those whose data states none
    info["not_perturbed"]       # held inputs the report names, as above
    info["has_gaps"]            # True if a nuclide, covariance block or spectrum
                                # had no usable uncertainty; of the held inputs
                                # above only the blocks in skipped_cross_material
                                # and skipped_nc set it
```

A ``sigma`` of ``0.0`` for a nuclide listed in ``no_covariance_data`` is not a
claim of certainty. It means the evaluation says nothing this run could use,
and the same goes for a half-life or decay energy listed as having no stated
sigma.

##### Uncertainty on activity, decay heat and dose

Activity, decay heat, contact dose and the decay-photon spectrum are all functions
of a *whole* inventory, so their spread has to be taken over the ensemble rather
than built from per-nuclide sigmas. Each has a getter that does exactly that,
evaluating the quantity once per replica and taking the spread of the result:

<!-- doctest: skip -->
```python
activity = results.get_activity_uncertainty(material_id=mid, step=1)   # Bq
heat = results.get_decay_heat_uncertainty(material_id=mid, step=1)     # W
dose = results.get_contact_dose_uncertainty(material_id=mid, step=1)   # Gy/h

if activity is not None:
    print(f"{activity.nominal:.4g} +/- {activity.std_dev:.2g} Bq "
          f"({activity.relative_std_dev:.1%}, {activity.replicas} replicas)")
```

Each returns an ``Estimate``. ``nominal`` is the unperturbed run and is there
whether or not uncertainty was asked for, while ``mean``, ``std_dev`` and
``relative_std_dev`` are ``None`` below two replicas, because a spread over fewer
than two samples is unmeasured rather than zero. Comparing ``mean`` against
``nominal`` is worth doing, with care over what a gap means. With the cross
sections alone and little burnup, the inventory is close to linear in the drawn
rates, whose mean is the nominal rate, so the two agree to within the sampling
error. With half-lives perturbed they need not: a nuclide's activity after a long
cooling is convex in its drawn half-life, so its mean sits above nominal although
the draw itself is unbiased, by nearly 20% for a 3% half-life sigma thirty
half-lives into cooling. Where neither applies, a gap says the perturbation is
biased rather than merely wide, and the floored counters above say where a draw
was cut.

``by_nuclide=True`` returns a ``dict[str, Estimate]`` instead of one total. Those
do not add up to the total in quadrature, and are not meant to.

Activity, decay heat and the photon spectrum need the material's ``volume``, since
a quantity in becquerel does. Contact dose does not, because the estimate takes
the material for a half-space.

Contact dose is also the one quantity here that is **not** linear in the atom
densities: the material attenuates its own decay photons, so the emitters sit in
the numerator and the whole inventory in the denominator, and a replica that makes
more of an emitter also absorbs more of it. Scaling the nominal dose by the
density spread would report the activity's uncertainty, which is wrong by the
whole value; evaluating per replica gets the cancellation for free.

The decay-photon spectrum comes back line by line:

<!-- doctest: skip -->
```python
lines = results.get_decay_photon_spectrum_uncertainty(material_id=mid, step=1)
for line in (lines or [])[:5]:
    print(f"{line.energy:10.1f} eV  {line.nominal:.3e} photons/s  "
          f"emitted in {line.emitting}/{line.replicas} replicas")
```

The set of lines is not the same in every replica: a nuclide that falls below the
density floor in one draw takes its lines out of that draw. The spectrum is
reported over the union, a line missing from a replica counts as a zero in it, and
``emitting`` says how many replicas emitted it at all, which is the difference
between a line that is dim and a line that is sometimes not there.

For any other quantity, ``get_uncertainty_inventories`` hands back every replica's
full inventory to evaluate over:

<!-- doctest: skip -->
```python
inventories = results.get_uncertainty_inventories(material_id=mid, step=1)
# evaluate the quantity once per replica, then take its spread
```

Summing per-nuclide sigmas in quadrature instead would assume the nuclides vary
independently, and they do not: a parent and its daughter move together, because
every daughter atom came out of a parent atom.

See also [Model.simulate_transmutation()](simulate_transmutation.md) for coupled transport-transmutation.
