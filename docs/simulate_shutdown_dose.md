# Shutdown Dose Rate

After irradiation, activated materials keep emitting decay gammas; the dose rate once the neutron source is switched off is the *shutdown dose rate*. YAMC computes it from a **single** transport run, one of the single-transport methods that also include UKAEA's **D1S** (Direct-1-Step) and **N1S** (Novel-1-Step), as opposed to **R2S** (Rigorous-2-Step), which runs a separate activation step. Decay photons are coupled directly into the neutron transport run, replacing the prompt photons. YAMC's method derives from D1S but computes time correction factors (TCFs) automatically for every reachable radionuclide and applies them in post-processing, so the irradiation and cooling schedule can change without re-running transport.

### Overview

1. Bin a photon tally by parent radionuclide (``parent_nuclides=model.radionuclides()``)
2. Run transport with ``use_decay_photons=True``
3. Describe the irradiation/cooling timeline as a ``PulseSchedule``
4. ``schedule.time_correct_tally(...)`` -> dose rate(s)

### Workflow

<!-- doctest: skip -->
```python
import yamc

_chain = "transmutation-endf-b8.1-sfr.arrow"
yamc.transmutation_decay_data = _chain
yamc.transmutation_reactions = _chain
yamc.transmutation_fission_yields = _chain

# 1. Bin a photon tally by parent radionuclide.
#    parent_nuclides=... is what makes a tally decay-photon-ready: it scores
#    one bin per radionuclide so the per-nuclide decay photons can be scaled
#    independently in post-processing. The radionuclides come straight from
#    the model's materials, so no schedule is needed yet. Works on any spatial
#    filter -- here a mesh, but cells/materials work too.
mesh = yamc.RegularRectangularMesh.from_domain(geometry, shape=1_000_000)
dose_tally = yamc.Tally(
    scores=["flux"],
    particle="photon",
    parent_nuclides=model.radionuclides(),
    mesh=mesh,
)

# 2. Run transport with decay photons
model = yamc.Model(
    geometry=geometry,
    tallies=[dose_tally],
    source=source,
    transport_secondary_photons=True,
    use_decay_photons=True,
)
results = model.simulate_transport(total_particles=1000000)

# 3. Describe the irradiation/cooling timeline. Each Pulse carries its own
#    source + rate (n/s) + (value, unit) duration; a Cooldown has no source.
schedule = yamc.PulseSchedule([
    yamc.Pulse(source=source, rate=1e14, duration=(1, "a")),   # 1 year irradiation
    yamc.Cooldown(duration=(30, "d")),                          # 30 days cooling
])

# 4. Time-correct the tally into dose rate(s). steps defaults to all schedule
#    steps; an int selects one (negative counts from the end). There is no
#    pre-irradiation baseline row to slice off.
dose = schedule.time_correct_tally(results[dose_tally], steps=-1)  # end of cooling
print(f"Dose rate at end of cooling: {dose.mean}")   # one value per mesh cell
print(dose.std_dev)                                  # matching uncertainty
print(dose.by_nuclide)                               # per-radionuclide breakdown
```

For a single ``int`` ``steps`` the ``DoseResult`` arrays are 1-D; for the
default (all steps) or a list of steps they are 2-D (one row per step), and
``dose.times`` gives the elapsed time at the end of each selected step. The
flat bin layout per row is ``score → mesh`` (the general order is in
[Spatial binning](tallies.md#spatial-binning)), so reshape with
``np.array(dose.mean).reshape(mesh.shape)``.

### Mixed-spectrum (DD + DT) campaigns

A schedule can mix sources -- e.g. a DT pulse, cooling, then a DD pulse. Each
distinct spectrum activates the materials differently, so run transport once
per source and pass ``time_correct_tally`` a dict mapping each source in
``schedule.sources`` to its run's tally result; the per-source corrections are
summed exactly.

<!-- doctest: skip -->
```python
dt = yamc.NeutronSource(energy=yamc.sources.fusion_neutron_spectrum(reactants="DT"))
dd = yamc.NeutronSource(energy=yamc.sources.fusion_neutron_spectrum(reactants="DD"))

schedule = yamc.PulseSchedule([
    yamc.Pulse(source=dt, rate=1e18, duration=(1, "a")),
    yamc.Cooldown(duration=(1, "d")),
    yamc.Pulse(source=dd, rate=5e17, duration=(30, "d")),
    yamc.Cooldown(duration=(7, "d")),
])

results_dt = yamc.Model(..., source=dt, ...).simulate_transport(total_particles=1_000_000)
results_dd = yamc.Model(..., source=dd, ...).simulate_transport(total_particles=1_000_000)
dose = schedule.time_correct_tally({dt: results_dt[dose_tally],
                                    dd: results_dd[dose_tally]})
```

### Reference

- ``model.radionuclides()`` -- activation radionuclides reachable from a Model's materials; pass to ``Tally(parent_nuclides=...)``.
- ``yamc.Pulse(rate, duration, source=None)`` -- an irradiation campaign. ``duration`` is seconds or a ``(value, unit)`` tuple (units: s, min, h, d, a). ``source`` is the NeutronSource for transport / SDR; omit it (sourceless) for ``Material.transmute``, where ``rate`` is a flux multiplier.
- ``yamc.Cooldown(duration)`` -- a cooling step (no source).
- ``yamc.PulseSchedule([steps...])`` -- the timeline; ``.sources`` lists its distinct sources, ``len(schedule)`` its step count.
- ``schedule.time_correct_tally(results, steps=None)`` -> ``DoseResult`` with ``.mean`` / ``.std_dev`` / ``.by_nuclide`` / ``.times``.

Parent-nuclide binning itself is just the ``parent_nuclides=`` argument on
``Tally`` (step 1) -- no separate setup call is needed.
