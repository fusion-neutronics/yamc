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

Density is always a float; its meaning is set by the optional ``units`` keyword. With ``units="atom/barn-cm"`` the value is the total atom density and the composition values are relative atom fractions (requires ``fraction_type="atom"``).

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

``material.density`` is always a float: the mass density for ``g/cm3``/``kg/m3``, or the total atom density for ``atom/barn-cm`` (and for transmuted materials).

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
    # transmutable=False is the default; composition stays fixed
)
```

The flag is not needed for ``Material.transmute()`` (standalone transmutation with a pre-computed flux). That method always evolves the material it is called on.

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

A material can also be taken directly from a curated collection instead of typing in its composition. ``yamc.materials.pnnl`` holds 410 materials from the [PNNL Compendium of Material Composition Data for Radiation Transport Modeling (PNNL-15870, Rev. 2)](https://www.pnnl.gov/main/publications/external/technical_reports/PNNL-15870Rev2.pdf): concretes, steels, tissues, air, soils, shielding, scintillators, fuels and more. The data is bundled in the wheel, so a lookup needs no network and cannot fail.

Materials are keyed by the report's own name, verbatim, so an entry can be quoted straight back to it.

```python
steel = yamc.materials.pnnl["Steel, Stainless 304"]
```

The result is an ordinary ``Material``, built fresh on every lookup, so it can be renamed and modified freely. Rename it when a mesh geometry expects a particular material name:

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

An unrecognised name raises ``KeyError`` and suggests near matches.

```python
yamc.materials.pnnl["Steel, Stainless 999"]
# KeyError: no material named 'Steel, Stainless 999' in the pnnl collection;
#           did you mean: Steel, Boron Stainless, Steel, HT9 Stainless, ...
```

#### Overriding fields at lookup

``temperature``, ``volume``, ``transmutable``, ``name`` and ``id`` are all settable on the returned material. ``density`` is not, because a ``Material``'s per-nuclide atom densities are fixed at construction. Use ``material()`` to build at a different density, for example a concrete at a different porosity.

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

Values are the report's own atom fractions, taken verbatim from its Elemental column, so no weight-to-atom conversion is applied on the way in.

Where an element is isotopically natural the element symbol is stored, and yamc expands it with its own natural abundances instead of the report's. Where a material is not natural (enriched lithium and boron, uranium, plutonium, heavy water, He-3) the individual nuclides are stored instead.

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

The example below uses a simple hand-specified 3-group spectrum so it can be read in isolation. In a real workflow the spectrum is usually obtained from a transport tally. See [Neutron and photon flux spectra](tallies.md#neutron-and-photon-flux-spectra) for how to score an energy-binned flux on a cell and feed ``results[tally].mean`` (scaled by the source strength) into ``transmute()``.

<!-- doctest: skip -->
```python
# The transmutation network is assembled from the per-subsection sources;
# point the three network-forming parts at the chain file.
yamc.transmutation_decay_data = "transmutation-endf-b8.1-sfr.arrow"
yamc.transmutation_reactions = "transmutation-endf-b8.1-sfr.arrow"
yamc.transmutation_fission_yields = "transmutation-endf-b8.1-sfr.arrow"

# The spectrum rides on each irradiation pulse: a NeutronSource whose energy is
# a Histogram over the energy groups. The Histogram normalizes the shape, so a
# multigroup flux from a tally can be passed straight in; see "Neutron and
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
- ``schedule``: a ``PulseSchedule`` of ``Pulse`` / ``Cooldown`` steps. Each irradiation ``Pulse`` carries a duration (seconds, or a ``(value, unit)`` tuple), a ``rate`` (the total flux magnitude [n/cm²/s]), and a ``source``, which is a ``NeutronSource`` whose energy is a ``Histogram`` giving the multigroup spectrum (its boundaries may be a named group structure such as ``"CCFE-709"``). The source's position/direction are ignored (no transport runs). ``Cooldown`` is decay-only. Different pulses may carry different spectra.
The transmutation network is assembled from the configured per-subsection sources (``yamc.transmutation_decay_data``, ``yamc.transmutation_reactions``, ``yamc.transmutation_fission_yields``, and optionally ``yamc.transmutation_branch_ratios``).

#### Clearance

A transmuted material can be checked against the regulatory clearance,
exemption and disposal limits directly, from yamc-core 0.14.0.
``clearance_index()`` divides each radionuclide's activity by its limit in a set
and sums the ratios; the material meets the limits when the sum is below the
set's threshold, normally 1. The sets cover the UK (EPR 2016, IRR 2017),
Germany (StrlSchV), the US (NRC 10 CFR 61, Fetter), the EU Basic Safety
Standards and IAEA GSR Part 3, and come from the
[radiological-material-clearance-finder](https://github.com/fusion-energy/radiological_material_clearance_finder)
crate, whose documentation describes each one and where its numbers come from.

The half-lives are the configured transmutation chain's, the same ones
``activity()`` and ``decay_heat()`` use. The composition is enough for the Bq/g
and Ci/m3 sets; only the total activity sets need ``material.volume``. Stable
isotopes have to stay in the material, since they are most of its mass and so
most of the Bq/g denominator.

<!-- doctest: skip -->
```python
if final is not None:
    result = final.clearance_index("UK_EPR16_out_of_scope")
    result.index        # float, sum of activity-to-limit ratios
    result.clearable    # bool, index below the set's threshold
    result.by_nuclide   # dict[str, float], each nuclide's ratio, largest first
    result.uncovered    # dict[str, float], activity no limit in the set covers
    print(result)       # a short report of the dominant nuclides

    results = final.clearance_indices()   # dict[str, ClearanceResult], every set
    routes = [name for name, r in results.items() if r.clearable]
```

``metal=True`` applies the NRC's activated metal rows, ``apply_default_limit``
and ``exclude_daughters`` switch off the catch-all limit and the secular
equilibrium credit that some tables define, and ``uncovered_fraction`` is the
number to check before trusting a comfortable index.

#### Activity, decay heat and decay photons

A transmuted material provides the usual decay quantities directly. Each takes
``per=None`` for the whole material, which needs ``volume``, or ``"cm3"`` or
``"g"`` for a density, and the totals also take ``by_nuclide=True``.

<!-- doctest: skip -->
```python
if final is not None:
    final.activity(per="g")                           # Bq/g
    final.decay_heat(per="cm3")                       # W/cm3
    final.decay_heat(per="cm3", component="gamma")    # the gamma part alone
    energies, rates = final.decay_photon_spectrum(per="cm3")   # lines, photons/s/cm3
    for continuum in final.decay_photon_continua(per="cm3"):
        print(continuum.nuclide, continuum.emission_rate)       # photons/s/cm3
```

``component`` splits decay heat into the recoverable-heat components the
decay data gives (``"beta"``, ``"gamma"`` or ``"alpha"``): the gamma heat is the
part that leaves a thin component, while the beta and alpha heat stay where they
are made. It raises if a nuclide making decay heat carries no split, rather than
understating the component by that nuclide's share.

``decay_photon_spectrum()`` gives the discrete lines. Part of some decay spectra
is given as a density over energy instead: the spontaneous-fission photons of an
actinide, or the whole emission of a nuclide far from stability whose lines were
never measured. ``decay_photon_continua()`` returns those, one
``PhotonContinuum`` per nuclide and continuum, each with its own ``energies``,
``rates`` per eV, ``interpolation`` law and ``emission_rate`` (the integral).
Their rates are per eV, so they cannot be added to the line rates; a material
with no continuum returns an empty list.

#### Repeated calls, and letting the data go

``transmute()`` loads the cross sections it needs into the material and keeps
them there, so calling it again on the same material does no file reading at all.
On a steel against ENDF/B-VIII.1 that is 1.8 s of a 2.3 s call. The composition is
not modified.

This helps when a material is transmuted more than once. It hurts in a sweep over
thousands of distinct compositions, where every material would hold a few hundred nuclides' data for as long as it lives.
``release_nuclear_data()`` drops it. The material stays usable, and the next call
that needs the data loads it again.

<!-- doctest: skip -->
```python
for composition in many_compositions:
    scan_mat = yamc.Material(composition=composition, density=7.9)
    scan_mat.transmute(schedule=schedule)
    scan_mat.release_nuclear_data()
```

#### Many materials at once

``yamc.transmute()`` is the plural of ``Material.transmute()``, for the cells of a
mesh from a transport run, a component broken into regions, or a sweep over
compositions. Each material gets exactly the answer its own ``transmute()``
would give, and the result is one ``TransmutationResults`` keyed by each
material's ``id``. The chain is read once, each nuclide's cross sections are
decoded once and shared, and a multigroup collapse runs once for every distinct
combination of spectrum, composition, temperature and shielding, before the
materials are solved in parallel.

<!-- doctest: skip -->
```python
cells = [
    yamc.Material(composition={"Fe": 1.0}, density=7.874, id=i + 1, volume=1.0)
    for i in range(3)
]
cell_fluxes = [1e14, 5e13, 1e13]   # each cell's flux magnitude [n/cm^2/s]
schedules = [
    yamc.PulseSchedule([
        yamc.Pulse(rate=flux, duration=(1, "y"), source=spectrum),
        yamc.Cooldown(duration=(1, "d")),
    ])
    for flux in cell_fluxes
]
many = yamc.transmute(cells, schedules)

many.get_final_material(material_id=2)
many.get_source_rates(material_id=2)   # that cell's flux magnitude per step
many.collapse_reuse                    # {"performed": 1, "requested": 3}
```

Give one schedule for every material or one per material, in the same order.
The schedules may differ in their spectra and rates but must share one
timeline, the same durations with irradiation on the same steps. The ids must
be distinct. ``data_uncertainty``, ``self_shielding_shape`` and
``self_shielding_chord`` apply to every material as ``Material.transmute()``
applies them to one. ``collapse_reuse`` reports the collapses actually run
against one per distinct spectrum per material: above, the three cells share
one spectrum and one composition, so they share one collapse whatever their
flux magnitudes.

#### Resonance self-shielding

The spectrum handed to ``transmute()`` is the flux outside the material. Inside
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
longest chord and shields more than any other shape of the same size. That makes it
an upper bound, not a safe default, so there is no default shape.

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

A ``strongest_factor`` of ``1.0`` means the run shielded nothing in practice.
That is not the same as not having tried.

A dilute run fills ``would_shield`` instead: the nuclides whose own resonances
are structured enough to have suppressed a reaction, each mapped to how strongly.
It is an indicator, not a correction and not a bound. It is computed from the one
reaction and that nuclide's own density, with no geometry in it, so it can sit
either side of the real factor: a lump thinner than a centimetre of chord shields
less, and a material whose other nuclides dominate the total at the resonance dips
the flux further. On the FNS tungsten foil it reads 0.634 for W186, where the
slowing-down correction on the same foil and spectrum gives 0.730. It warns that
the answer may be high and that the nuclide is a resonance absorber. Read the size
of the effect off a shielded run.

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
material's own spectrum, and their resonance parameters where the library
carries them; the half-lives, decay energies, decay photon lines and decay
branching ratios from the decay data's own standard deviations; the fission
yields from theirs; the isomeric splits from MF=40; the flux from the error you
give it; and, on
[Model.simulate_transmutation()](simulate_transmutation.md), the
transport-tallied rates from their per-history covariance. The schedule is
re-solved with the drawn rates, branching, flux and half-lives, and the drawn decay energies
are applied where decay heat is evaluated. That is exact to all orders in the
matrix exponential: nothing is linearized, and the solver itself is untouched,
only its input changes.

``DataUncertainty(seed=...)`` is reproducible. A nuclide's cross-section,
half-life, decay-branching, decay-energy, decay-photon, fission-yield and
isomeric-branching draws in a given sample are a pure function of
``(seed, sample, nuclide)``, and the flux draw of ``(seed, sample, spectrum)``,
the spectrum's place among the schedule's spectra. The statistical draw is keyed
on ``(seed, sample)`` and made jointly over the material's tallied rates, so
adding or removing a tallied rate changes the draws of the others. It carries no
material key, so every material in one run draws the same deviates in a given
sample: each material's own sigma is right, but the statistical draws of
different materials are not independent, and replica inventories should not be
combined across materials. The flux draw has the same limit: it is keyed on the
spectrum's place in each material's own schedule, not on the material, so on a
multi-material ``transmute()`` the materials whose schedules put a spectrum in the
same place draw identical flux deviates. Their flux errors are then fully
correlated, which is right for one shared tally and wrong for separately tallied
cells. In every case the same seed gives the same answer whatever the sample
count or iteration order. Leave ``samples`` unset and the solver adds samples
until the reported sigmas settle; setting it bounds cost or reproduces a
specific run, but it is not an accuracy dial.

Omitting the argument costs nothing. No covariance is read from disk, no matrix
is folded or factorized, and the inventories are bit-identical.

##### Choosing what to perturb

Nine sources can be perturbed, and by default every one that applies to the call
is:

<!-- doctest: skip -->
```python
yamc.DataUncertainty.available_sources()
# ['cross_sections', 'flux_spectrum', 'half_life', 'decay_branching',
#  'statistical', 'decay_energy', 'decay_photon_lines', 'fission_yield',
#  'isomer_branching']

# isolate one contribution
xs_only = yamc.DataUncertainty(seed=42, sources=["cross_sections"])
```

- ``cross_sections``: the activation cross sections, from the evaluation's MF=33
  covariance folded against the spectrum the rates are collapsed with. On
  ``Model.simulate_transmutation`` the fold uses the tallied spectrum,
  normalised the same way as the rates it divides by. Where a nuclide's data
  folder carries its resonance parameters (MF=2 and MF=32,
  ``resonance_parameters.arrow``, in the published libraries), each sample
  draws the parameters themselves, rebuilds the resonance range's cross
  sections from the draw and Doppler broadens the change, so that range's
  uncertainty is exact in the parameters rather than first order. On
  ``Model.simulate_transmutation`` the transport also tallies the flux on each
  such nuclide's own energy grid through its resonances, and the change is
  folded against that, self-shielding dips included. A nuclide whose
  parameters cannot be sampled keeps the first-order MF=32 rows of
  ``covariance.arrow``, and the report's ``resonance_parameters`` says which
  and why.
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
  density's whole spread. A shutdown-dose ``PulseSchedule.time_correct_tally``
  applies it, and ``decay_photon_lines`` (below), and no other source.
- ``decay_branching``: the decay branching ratios of each reachable parent with
  two or more modes and a stated sigma, however many modes it has. Each draw
  keeps the parent's total and every ratio between zero and it, and each mode
  takes its evaluated ratio as its mean and its stated sigma as its standard
  deviation, so a small branch with a large relative sigma is drawn skewed
  rather than held. A mode is held only when its stated sigma is one no ratio
  with that mean can have, or it states none and is not the complement of the
  rest. The report names every parent held whole
  (``no_decay_branching_uncertainty``, ``decay_branchings_too_wide``) or drawn
  short of its stated sigmas (``decay_branchings_partly_held``,
  ``decay_branchings_sigma_adjusted``). Only the inventory moves: a drawn
  parent's lines and decay energy per decay follow its nominal branching.
- ``statistical``: the Monte Carlo uncertainty of transport-tallied reaction
  rates, from their per-history covariance, so rates scored by the same
  histories move together.
  [Model.simulate_transmutation()](simulate_transmutation.md) only, which
  accepts ``data_uncertainty`` with ``method="independent"`` alone:
  ``method="coupled"`` raises.
- ``decay_energy``: each nuclide's mean decay energy, from the sigmas the decay
  data gives its recoverable-heat components (beta, gamma, alpha), or from the
  total's sigma where no component carries one. The evaluation
  states no correlation between the components, so they are drawn as
  independent and the total's sigma is taken as their quadrature sum. A decay
  energy never enters the solve, so this moves decay heat and nothing else.
- ``decay_photon_lines``: each decay photon spectrum's normalisation (FD for
  lines, FC for a continuum), one draw per spectrum shared by all its lines, and
  each line's own intensity (dRI) and energy (dER), from the decay data's
  MT=457 sigmas. It moves the decay photon spectrum, the contact dose and,
  through the gamma decay energy, the decay heat; no photon enters the solve.
  The decay data does not state how a nuclide's photon intensities correlate
  beyond the shared normalisation, so each photon output is given at both
  ends: ``std_dev`` takes them independent and ``std_dev_correlated`` fully
  correlated.
- ``fission_yield``: each fissioning parent's independent yields (MT=454), from
  the sigma the evaluation states on each, drawn on the tape's own products and
  summed onto the chain's the way the converter summed the nominal yields. No
  evaluation states a correlation between yields, so each is drawn
  independently and a draw's yields are not renormalised. It needs a chain
  carrying the evaluated yields (``fission_yields/evaluated_yields.arrow``);
  without them every parent is listed under ``no_fission_yield_uncertainty``.
- ``isomer_branching``: the isomeric split of every reachable reaction whose
  evaluation states the MF=40 covariance of its final states' production (the
  branching subsection's ``branching_covariance.arrow``). Each state's
  production is drawn on the covariance's own grid, folded against the same
  spectrum as the nominal split, and the split remade from the draw. Where the
  list names its ground state only the shares move, and the reaction total
  stays the cross sections' own, so the total's uncertainty is not counted
  twice. An isomer listed alone, and every ``(n,n')`` isomer, is made at its
  drawn production, so its MF=40 sigma replaces the MF=33 sigma of its reaction
  that it carries under ``cross_sections`` alone, and can be the smaller. The
  spectrum calls only: a transport run's split comes from tallied partials and
  is held.

A call ignores the sources that do not apply to it, and ``sources`` in the report
below lists the ones that did.

Restricting ``sources`` isolates one contribution, so the growth of the inventory
sigma as a source is added can be measured directly. Naming a source this build
cannot perturb raises an error instead of being ignored, so a source that is not
implemented yet cannot look like one that contributed nothing.

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

##### Where the uncertainty comes from

``DataUncertainty(attribution=True)`` also says which source and which
evaluation a nuclide's spread comes from, read with
``get_uncertainty_breakdown``. It is off by default because it costs further
solves, one ensemble per source and one deterministic solve per contributor,
and it changes none of the numbers the run otherwise reports.

<!-- doctest: skip -->
```python
attributed = material.transmute(
    schedule=schedule,
    data_uncertainty=yamc.DataUncertainty(seed=42, attribution=True),
)
breakdown = attributed.get_uncertainty_breakdown(material_id=mid, nuclide="Mn56", step=1)
if breakdown is not None:          # None unless attribution=True
    breakdown["variance"]          # the total, the square of get_nuclide_uncertainty
    breakdown["by_source"]         # each source alone, resampled; these sum to the total
    breakdown["unattributed"]      # what that sum leaves: interaction and sampling noise
    breakdown["contributors"]      # [(source, nuclide, reaction, variance)], largest first
```

``by_source`` resamples each source on its own, so it says how much of the
spread is statistical and how much is each kind of nuclear data; the sources
are independent and their variances add. ``contributors`` is first order and
names the evaluation to look at: within the cross sections a nuclide's whole
evaluation has a ``reaction`` of ``None`` and each channel alone names its
reaction, and a half-life or a decay branching has ``None`` too. The
total stays the resampled one.

##### How a sample is drawn

A cross-section channel's rate is multiplied by a lognormal factor with mean one
and that channel's variance in the matrix it samples from (the folded covariance,
or the repaired one described below), driven by the channel's correlated Gaussian
deviate. So each channel keeps that mean and variance exactly and cannot go
negative. Between channels the draw keeps the ordering of the correlated Gaussian
deviates, so their rank correlation, but not the stated linear correlation,
which shifts by an amount that vanishes as the sigmas go to zero and grows with
them. For a channel whose relative sigma approaches or exceeds 1 the lognormal
is heavy-tailed, so an ensemble of the size the solver stops at (at
most 1024 unless ``samples`` is set) tends to report a sigma below the stated
one, and the settling test cannot detect that.

A half-life or decay energy is multiplied by a lognormal factor of the same
kind, with mean one and variance equal to its squared relative sigma, drawn per
nuclide (for decay energy, per beta, gamma and alpha component, or on the total
where no component carries a sigma). The decay data states an expected value and
a standard deviation and no correlation, so each draw keeps exactly what the
evaluation states and cannot go negative, and there is no floor. A shutdown-dose
``time_correct_tally`` draws its half-lives the same way. A sigma no draw can
carry, one stated on a decay energy of zero or one that is not finite, is held at
nominal and listed in ``half_life_uncertainty_not_carried`` or
``decay_energy_uncertainty_not_carried``.

The flux and tallied-rate sources draw from a Gaussian on the stated sigma or
covariance, since they carry correlations a lognormal would not keep. A flux bin
or tallied rate drawn below zero is set to zero, and
``get_data_uncertainty_info`` counts those (``flux_bins_floored`` and
``statistical_floored``, each beside its ``_sampled`` total). A count well above
zero says the sigma is a large enough share of the value that a Gaussian no
longer describes it.

A parent's decay branching ratios sum to a fixed total, so the deviations of a
draw sum to zero. Of the zero-sum covariances that realise every stated sigma,
the draw takes the one of largest entropy, which assumes no correlation the sum
rule does not force (for two or three modes it is the only one); where the
largest sigma exceeds the sum of the others', that mode is drawn at their sum.
Every moving mode but one is a logit-normal on ``[0, T]``, ``T`` the parent's
total, fitted so that its mean and standard deviation are the evaluated ratio
and its sigma, and the modes are coupled through a Gaussian copula that carries
the target covariance. The last mode is their complement, which keeps the
total. A draw whose complement would fall below zero is pulled back towards
nominal and counted in ``decay_branchings_floored``, beside
``decay_branchings_sampled``.

Evaluated covariance matrices are not always positive semi-definite, and one that
is not has no square root to sample from. The negative eigenvalues of the folded
matrix are set to zero, which is the standard repair, and the repair adds
variance rather than removing it: no channel's sampled sigma falls below the
evaluated one, and a channel the evaluation gave no variance can pick some up.
A matrix counts as repaired when its correlation matrix, over the channels with
a positive stated variance, has an eigenvalue below ``-n * 1e-12`` (``n`` the
number of those channels), or when a channel is stated with a negative variance,
or with a zero one and a covariance to another channel. Judged on the
correlation matrix, the round-off of the eigendecomposition is on each channel's
own scale, so a PSD matrix that comes back a few ulps negative is not a repair,
and a channel with a huge variance cannot hide a real repair among the others.
``covariance_repaired`` lists the distinct nuclides, not a count over spectra,
with a repaired channel a draw can move (a positive rate on a spectrum the
schedule irradiates with), and any makes ``has_gaps`` true. The report covers
the nuclides the material can populate: the fold takes every chain nuclide with
data, which from almost any composition is the chain's whole closure, so a
nuclide is reported only when an upper bound on its density over the schedule,
at nominal rates, reaches the solver's floor of ``1e-30`` atoms/b-cm. One the
bound leaves out cannot move any nominal density by as much as that floor, but a
replica's rates on a wide channel can sit orders above nominal, so the repaired
nuclides outside the bound with a channel a draw can move are named in
``covariance_repaired_outside_bound``, and any of those makes ``has_gaps`` true
as well, since the bound says nothing about a replica. ``covariance_repairs``
gives one entry per repaired populated nuclide and spectrum, including repairs
no draw can move, with ``lambda_min``, ``lambda_max``, ``clipped_fraction``
(the variance the repair added over the stated trace) and, per channel,
``evaluated_variance`` (the folded diagonal as stated, which can be negative),
``evaluated_sigma`` (``None`` for a negative variance) and ``sampled_sigma``.
A small ``|lambda_min| / lambda_max`` does not mean a small widening of the
channel that matters, which is why each channel's two sigmas are given.
``worst_sigma_inflation`` is the largest sampled over evaluated sigma, minus
one, over the repaired channels a draw can move, and
``rate_weighted_sigma_inflation`` is the weighted mean of each channel's own
sampled over evaluated sigma, minus one, over every sampled channel of a
populated nuclide, each weighted by its rate at unit flux times its spectrum's
fluence in the schedule times its parent's initial density. The weight is on
each channel's inflation, not on its sigma, so a wide channel with a small rate
cannot drown out a repair on the channels that carry the reactions. The initial
density makes it a measure of first-generation reactions: a nuclide the material
starts without carries no weight, and its repairs show in
``worst_sigma_inflation`` and ``covariance_repairs`` instead. The sampled sigma
is read off the factorization for every matrix, so on a matrix that needed no
repair the decomposition's round-off shows here as it is, far below any real repair.
Either headline can be ``float('inf')`` when a repair gives spread to a channel
stated with none, and strict JSON does not accept that.

``sigma_at_least_one`` and ``sigma_at_least_ten`` name the sampled channels,
keyed ``"Nuclide (n,x)"``, of populated nuclides with a positive rate on a
spectrum the schedule irradiates with, whose folded relative sigma as evaluated
is at least one or ten. At that width the spread depends on the lognormal chosen to carry
the evaluation's two moments, not on the evaluation alone.
``sigma_at_least_one_outside_bound`` does the same for the nuclides outside the
populated bound, since a replica's draw on exactly such a channel can sit orders
above nominal and populate the nuclide; the ten-or-more subset reads off its
values. They are not a gap.

##### What it does not cover

Every input below is held at its evaluated or nominal value in every sample, or
drawn without a correlation the data states. Some carry published uncertainties
the converted data does not hold yet, some have none to read, and some are a
property of the model rather than of the data. Most leave a term out of the
reported spread, so the sigma is too small by it. Some can go either way:
self-shielding and the flux's response to the data, because holding them also
removes a feedback that narrows the spread, and the dropped correlations, because
a missing correlation can raise or lower the sigma of a sum. Each item says
which way it goes, and where.

- **Decay branching ratios the draw holds.** A mode whose stated sigma no
  ratio with its mean can have, a mode stating no sigma that is not the
  complement of the rest, and a parent with no sigma at all are held at their
  evaluated ratios, and each is listed in the report under the reason. A drawn
  parent's per-decay emission, its lines and decay energy, stays that of the
  nominal scheme. The sigma is read from the decay data's
  ``branching_ratio_uncertainty`` column, so a library converted before that
  column existed carries none and lists its parents under
  ``no_decay_branching_uncertainty``.
- **Fission yield correlations.** The ``fission_yield`` source draws each
  stated yield independently, since no evaluation publishes a correlation
  between yields, and does not renormalise a draw's yields to the nominal
  total. A parent whose tape yields do not sum back to the chain's
  (``fission_yields_mapping_mismatch``) is held. The weights that combine a
  nuclide's yield sets by incident energy come from the nominal spectrum in
  every sample.
- **Isomeric branching without MF=40.** The ``isomer_branching`` source draws
  the split between a product's ground and metastable states where the
  evaluation states MF=40 on the spectrum calls. A split with no MF=40 on any
  of its states is held and listed in ``no_isomer_branching_uncertainty``, every
  split is held on a transport run, and MF=40 blocks that cannot be used are
  counted per parent in ``isomer_branching_blocks_not_used``. ENDF/B-VIII.1,
  FENDL-3.2d and TENDL state no correlation in MF=40 between a product's states,
  or with the reaction's MF=33, so the split's uncertainty is only partly
  stated ([Covariance is presence, not
  coverage](nuclear_data_libraries.md#covariance-is-presence-not-coverage)). With
  ``transmutation_branch_ratios`` configured, the ``(n,n')`` channels the
  overlay adds to reach a metastable state, In115 to In115m among them, come
  from its MF=10 partials with no group-averaged cross section behind them. They
  draw no cross-section covariance, and on the spectrum calls no
  ``flux_spectrum`` perturbation either, so they are missing from
  ``rate_fraction_covered``.
- **Resonance-parameter covariance (MF=32) neither sampled nor written.** Many
  evaluations keep the resolved-resonance uncertainty on capture, elastic and
  fission in MF=32 and leave MF=33 at zero over that range, which reaches down
  to thermal. Where the data folder carries the resonance parameters, the
  ``cross_sections`` source samples them exactly; otherwise the converter's
  first-order MF=32 rows in ``covariance.arrow`` are folded. What is left held
  is a range whose formalism is not reconstructed, and data converted before
  either was written, on which a product made by capture at thermal or
  resonance energies draws little or none of its cross section's uncertainty.
  ENDF/B-VIII.1 W186 is the example: its resonance parameters carry about 1.5%
  on capture below 10 keV, where its MF=33 gives zero, so with MF=33 alone the
  nuclide is still listed under ``perturbed`` and the gap shows in that
  channel's ``rate_fraction_covered``, which counts only rate from energies
  with a nonzero stated variance: W186 ``(n,gamma)`` reads about 0.07 on a
  fusion spectrum. On ``Model.simulate_transmutation`` the sampled parameters
  are folded against the nominal transport's flux, so for a self-shielded
  resonance the sigma is an overestimate: a stronger resonance would deepen its
  own dip. ``simulate_transport(data_uncertainty=...)`` carries that response
  through its replica weights.
- **MF=33 blocks that are not used.** What is folded is each reaction's own
  MF=33 blocks, including those whose ``MAT1`` names the evaluation itself, the
  blocks of a lumped reaction (MT=851 to 870) with one component, which are that
  component's, and the covariance an LTY=0 NC block derives from the reactions
  it names (ENDF/B-VIII.1 O16 ``(n,p)`` is stated only that way). Four kinds
  are left out. Blocks that correlate with another evaluation (``MAT1`` naming
  another material, the links to the standards among them) are skipped, so each
  nuclide's cross sections are drawn independently of every other nuclide's,
  and so are blocks whose partner is not a cross section (``XMF1`` other than 0
  or 3). NC blocks that cannot be derived (LTY 1 to 4, or an LTY=0 block that
  fails, counted in ``skipped_nc``) are skipped, so a reaction whose covariance
  is given only that way draws none. A lumped reaction with several components
  states the uncertainty of their sum and of none of them, so unless a
  derivation names it, it is not given to any component and is listed in
  ``lumped_covariance_not_assignable``. Tungsten's ``(n,2n)``, the channel
  behind most of a tungsten first wall's decay heat, is one: ENDF/B-VIII.1,
  JEFF-4.0 and FENDL-3.2d state its covariance only lumped with ``(n,2np)``.
  Blocks on partial levels (MT=600 to 849 and 875 to 891) are folded when an
  LTY=0 NC block names them, and are otherwise neither used nor counted, since
  the chain drives the total and has no rate for a level alone.
  ``skipped_cross_material`` and ``skipped_other_file`` count the first kind
  per nuclide, for reactions the fold reaches. A nuclide whose only blocks are
  of these kinds is listed under ``no_covariance_data``, and a channel missing
  its covariance shows in ``rate_fraction_covered`` and a lower
  ``rate_fraction_covered_total``.
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
  Those figures fold the shielded partial rates, as the covariance fold does:
  on a self-shielded run the partials come from the same shielded flux as the
  rate they are divided by.
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
- **Decay photon spectrum covariance and continuum shape.** The
  ``decay_photon_lines`` source draws each spectrum's normalisation and each
  line's intensity and energy. MT=457 can also state a covariance between a
  spectrum's lines (LCOV), and no library ships one, so the lines of a spectrum
  are drawn independently apart from their shared normalisation, and the
  photon outputs are given at both ends of that range (``std_dev`` and
  ``std_dev_correlated``). A continuum's shape is held, and only its
  normalisation is drawn. A shutdown dose draws the normalisations alone.
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
  stated range, which includes commercially depleted lithium, spans a factor of
  four). The bundled table carries each isotope's stated uncertainty, or for
  the twelve elements the tables give as an interval, that interval, readable
  through ``yamc.data.natural_abundance_records()``. For an interval element the
  tables state no distribution within it. Every sample uses the lookup value,
  so none of that spread is propagated: the tables' representative value, or for
  an interval element a single best measurement, which for Li6 (0.07589) sits
  near the top of its interval, [0.019, 0.078].
- **Atomic masses.** Each nuclide's mass, from AME2020, converts between mass and
  atoms the same way in every sample. AME2020 states an uncertainty on every
  mass, and the bundled table does not carry it.
- **Element fractions and density.** These are your input and start every sample
  at the values you gave. There is no argument to give their uncertainty, so it
  is not perturbed, and propagating it means sampling compositions yourself, one
  run per composition. For dose from a trace impurity it can outweigh everything
  above: in a 316L-type steel with 0.1 wt% cobalt, irradiated for five years at
  1e14 n/cm^2/s in an illustrative first-wall spectrum, the contact dose at ten
  years' cooling is almost all Co60, about 70% of it made from the cobalt itself,
  so a 30% uncertainty on the cobalt content is about 21% on the dose.

The report's ``not_perturbed`` names what the run held at nominal. Every run
names the isomeric branching of reactions with no stated MF=40, the MF=33
blocks that are not used (cross-material, a partner that is not a cross section,
NC blocks that cannot be derived, and lumps of several components), MF=32
neither sampled nor in ``covariance.arrow``, the decay photon spectrum
covariance and continuum shape, the three dose tables and the build-up factor,
material composition and density, natural isotopic abundance and atomic mass
(AME2020). The rest depend on the run. Half-life, decay branching, decay
energy, decay photon lines, fission yield, isomeric branching and the
activation cross sections are named when ``sources`` leaves them out, and when decay branching is drawn for any parent, the per-branch decay
emission (lines and decay energy following the nominal branching) is named
instead. On a transport call the flux's response to perturbed cross sections is
named whenever the ``cross_sections`` source is on, since the one transport is
not rerun per sample. On a spectrum call the flux spectrum is named when no
spectrum's flux is drawn, and named as held for spectra without a sigma only
when some are drawn and some are not. On a transport call the tallied-rate
statistics are named when the ``statistical`` source is off or there are no
tallied rates to draw. The self-shielding correction is named whenever
self-shielding is set.

Because a reported sigma of zero can mean "well known" or "nothing published",
the report, ``get_data_uncertainty_info``, keeps a nuclide with no usable MF=33
(``no_covariance_data``: none in the evaluation, or none this run could use,
such as blocks all skipped as above or all for reactions the run does not drive)
apart from one that states it. A channel whose stated covariance is zero over
part of its range, as in the MF=32 case above, stays under ``perturbed``, and
``rate_fraction_covered`` gives the share of its rate from energies with a
nonzero stated variance:

<!-- doctest: skip -->
```python
info = results.get_data_uncertainty_info(material_id=mid)
if info is not None:            # None unless data_uncertainty was passed
    info["sources"]             # the sources that applied to this run
    info["perturbed"]           # nuclides that had usable MF=33 covariance
    info["no_covariance_data"]  # nuclides with none usable, absent or skipped
    info["rate_fraction_covered"]  # per channel, the share of the rate the fold
                                   # divides by (shielded on a self-shielded run)
                                   # from energies with a nonzero stated variance
    info["rate_fraction_covered_total"]  # ... over the run, by production; None on
                                         # a transport run or when a partials map
                                         # lists a channel
    info["skipped_cross_material"]    # blocks naming another evaluation, not used
    info["skipped_other_file"]        # blocks whose partner is not a cross section
    info["skipped_nc"]                # NC blocks that could not be derived
    info["lumped_covariance_not_assignable"]  # lumps of several components
    info["partials_above_rate"]  # channels whose partial rates add up to more
                                 # than the rate: sigma overstated
    info["partials_below_rate"]  # ... add up to less: sigma understated
    info["sigma_at_least_one"]   # channels with a relative sigma of 1 or more
    info["half_lives_perturbed"]      # nuclides whose half-life was drawn
    info["no_half_life_uncertainty"]  # ... and those whose data states no sigma
    info["decay_branchings_perturbed"]      # multi-mode parents whose split was drawn
    info["no_decay_branching_uncertainty"]  # multi-mode parents with no sigma
    info["decay_branchings_partly_held"]    # drawn with some modes held
    info["decay_branchings_sigma_adjusted"] # drawn short of a stated sigma
    info["decay_photon_lines_perturbed"]    # nuclides whose photon data was drawn
    info["fission_yields_perturbed"]        # fissioning parents whose yields were drawn
    info["isomer_branchings_perturbed"]     # reactions whose isomeric split was drawn
    info["no_isomer_branching_uncertainty"] # splits held, no MF=40 on any state
    info["resonance_parameters"]  # per nuclide: parameters sampled, or
                                  # first-order rows and why
    info["decay_energies_perturbed"]  # nuclides whose decay energy was drawn
    info["no_decay_energy_uncertainty"]  # ... and those whose data states none
    info["half_life_uncertainty_not_carried"]     # a stated sigma no draw can
    info["decay_energy_uncertainty_not_carried"]  # carry, held at nominal
    info["not_perturbed"]       # held inputs this run names, as above
    info["has_gaps"]            # True if a nuclide, covariance block or spectrum
                                # had no usable uncertainty, a block was skipped,
                                # malformed or of an unsupported layout, either
                                # partials check listed a channel, a covariance
                                # was repaired, a half-life or decay energy had
                                # no stated sigma or one no draw can carry, a
                                # parent's decay branching was held or drawn
                                # short, or an isomeric split had no MF=40
```

A ``sigma`` of ``0.0`` for a nuclide listed in ``no_covariance_data`` is not a
claim of certainty. It means the evaluation says nothing this run could use,
and the same goes for a half-life or decay energy listed as having no stated
sigma.

##### Uncertainty on activity, decay heat and dose

Activity, decay heat, contact dose and the decay-photon spectrum are all functions
of a whole inventory, so their spread has to be taken over the ensemble rather
than built from per-nuclide sigmas. Each has a getter that does this,
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
``relative_std_dev`` are ``None`` below two replicas, because a spread over
fewer than two samples is unmeasured rather than zero. Comparing ``mean``
against ``nominal`` is useful, but a gap between them needs interpreting. With the
cross sections alone and little burnup, the inventory is close to linear in the
drawn rates, whose mean is the nominal rate, so the two agree to within the
sampling error. With half-lives perturbed they need not: a nuclide's activity
after a long cooling is convex in its drawn half-life, so its mean sits above
nominal although the draw itself is unbiased, by nearly 20% for a 3% half-life
sigma thirty half-lives into cooling. A product reached through two or more
perturbed steps that share a draw (the flux or the statistical draw) is
nonlinear in that draw too: its mean moves from nominal by about the covariance
of the shared draws, the product of the steps' relative sigmas times their
correlation, with no bias in the draw. That is above nominal for positively
correlated steps, such as two rates folded from one flux given by per-bin
sigmas, and below it where a stated correlation is negative. Outside these
cases, a gap says the perturbation is biased rather than merely wide, and the
floored counters above say where a draw was cut.

``by_nuclide=True`` returns a ``dict[str, Estimate]`` instead of one total. Those
do not add up to the total in quadrature, and are not meant to.

Activity, decay heat and the photon spectrum need the material's ``volume``, since
a quantity in becquerel does. Contact dose does not, because the estimate takes
the material for a half-space.

Contact dose is also the one quantity here that is not linear in the atom
densities: the material attenuates its own decay photons, so the emitters sit in
the numerator and the whole inventory in the denominator, and a replica that makes
more of an emitter also absorbs more of it. Scaling the nominal dose by the
density spread would report the activity's uncertainty, which is wrong by the
whole value; evaluating per replica captures the cancellation.

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
``emitting`` says how many replicas emitted it at all, which separates a dim line
from one that is sometimes absent.

For any other quantity, ``get_uncertainty_inventories`` returns every replica's
full inventory to evaluate over:

<!-- doctest: skip -->
```python
inventories = results.get_uncertainty_inventories(material_id=mid, step=1)
# evaluate the quantity once per replica, then take its spread
```

The inventories are densities only: a replica's drawn half-lives and decay
energies are not returned with them. Evaluating a quantity that reads a decay
constant or a decay energy from them with the chain's nominal values gives it
the wrong spread, the one the ``half_life`` item under [Choosing what to
perturb](#choosing-what-to-perturb) describes: a saturated activity then picks
up the density's whole spread. Activity, decay heat, contact dose and the photon
spectrum have the getters above. For any other such quantity, run with
``sources`` leaving out ``half_life`` and ``decay_energy``, and read its spread
as excluding those two.

Summing per-nuclide sigmas in quadrature instead would assume the nuclides vary
independently, and they do not: a parent and its daughter move together, because
every daughter atom came out of a parent atom.

See also [Model.simulate_transmutation()](simulate_transmutation.md) for coupled transport-transmutation.
