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
On a steel against ENDF/B-8.1 that is 1.8 s of a 2.3 s call. The composition is
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

``self_shielding_info`` reports what was done, and is always present:

<!-- doctest: skip -->
```python
shielding = results.self_shielding_info
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
could have suppressed a reaction, each mapped to the strongest suppression it
could have seen. That is a bound computed from the one reaction alone, ignoring
the rest of the material and the geometry, both of which push the real factor back
toward one. So it says "this answer may be high, and here is by how much at the
very most", which is the warning a dilute run of a resonance absorber should carry
rather than silence.

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

The activation cross sections are resampled from their ENDF MF=33 covariance,
folded against this material's own spectrum, and the schedule is re-solved for
each sample. That is exact to all orders in the matrix exponential: nothing is
linearized, and the solver itself is untouched -- only its input changes.

``DataUncertainty(seed=...)`` is reproducible: a nuclide's perturbation in a
given sample is a pure function of ``(seed, sample, nuclide)``, so the same seed
gives the same answer whatever the sample count or iteration order. Leave
``samples`` unset and the solver adds samples until the reported sigmas settle;
setting it bounds cost or reproduces a specific run, but it is not an accuracy
dial.

Omitting the argument costs nothing. No covariance is read from disk, no matrix
is folded or factorized, and the inventories are bit-identical.

##### Choosing what to perturb

Two sources can be perturbed, and by default both are:

<!-- doctest: skip -->
```python
yamc.DataUncertainty.available_sources()   # ['cross_sections', 'flux_spectrum']

# isolate one contribution
xs_only = yamc.DataUncertainty(seed=42, sources=["cross_sections"])
```

Restricting ``sources`` is how a run isolates one contribution, so that adding a
source and watching the inventory sigma grow is a measurement rather than a guess.
Naming a source this build cannot perturb raises rather than being ignored: a
source that has not landed yet must not look like one that contributed nothing.

``flux_spectrum`` propagates the error bar on the spectrum you supplied, which
means the pulse has to carry one. ``Pulse(flux_std_dev=...)`` takes the per-bin
standard deviation, in the same units and the same order as the ``Histogram``
values; from a yamc run that is ``results[tally].standard_deviation``.

<!-- doctest: skip -->
```python
pulse = yamc.Pulse(
    rate=total_flux,
    duration=(1, "h"),
    source=spectrum,
    flux_std_dev=flux_sigma,   # per bin, same order as the Histogram values
)
```

Omitting it is the common case: a spectrum taken from a published reference set
carries no stated error. The ``flux_spectrum`` source then contributes nothing and
says so in ``data_uncertainty_info["spectra_without_flux_sigma"]``, rather than
reading as a flux known exactly.

##### What it does not cover

Half-lives, decay branching ratios, fission yields and the isomeric-branching
overlay are held at their evaluated values, and they carry uncertainties of their
own that this does not propagate. ``data_uncertainty_info["not_perturbed"]`` names
them for the run in hand.

Because a reported sigma of zero can mean "well known" or "nothing published",
the two are kept apart in ``data_uncertainty_info``:

<!-- doctest: skip -->
```python
info = results.data_uncertainty_info
if info is not None:            # None unless data_uncertainty was passed
    info["perturbed"]           # nuclides that had usable MF=33 covariance
    info["no_covariance_data"]  # nuclides whose evaluation carries none
    info["rate_fraction_covered"]  # per channel, the share of the rate the
                                   # covariance grid actually spans
    info["not_perturbed"]       # the sources this does not propagate at all
    info["has_gaps"]            # True if anything above was left out
```

A ``sigma`` of ``0.0`` for a nuclide listed in ``no_covariance_data`` is not a
claim of certainty. It means the evaluation says nothing.

The report also carries ``matrices_clipped`` and ``worst_relative_clip``:
evaluated covariance matrices are not always positive semi-definite, and one
that is not has no square root to sample from. The negative eigenvalues are set
to zero, which is the standard repair, and how large a repair it was is recorded
rather than absorbed. ``rates_floored`` counts samples whose cross section went
negative and was truncated at zero; a large share means the Gaussian is being
used past the point where it describes the quantity.

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
``nominal`` is worth doing: activity and decay heat are linear in the atom
densities, so the two agree to within the sampling error, and a gap between them
says the perturbation is biased rather than merely wide.

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
