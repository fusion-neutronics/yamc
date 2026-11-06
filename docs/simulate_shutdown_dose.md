# Shutdown Dose Rate

After irradiation, activated materials keep emitting decay gammas. The dose rate once the neutron source is switched off is the shutdown dose rate. YAMC computes it from a single transport run. This puts it among the single-transport methods, which also include UKAEA's D1S (Direct-1-Step) and N1S (Novel-1-Step), as opposed to R2S (Rigorous-2-Step), which runs a separate activation step. Decay photons are coupled directly into the neutron transport run, replacing the prompt photons. YAMC's method derives from D1S but computes time correction factors (TCFs) automatically for every reachable radionuclide and applies them in post-processing, so the irradiation and cooling schedule can change without re-running transport.

### Overview

1. Bin a photon tally by parent radionuclide (``parent_nuclides=`` the model's ``radionuclides()``)
2. Run transport with ``use_decay_photons=True``
3. Describe the irradiation/cooling timeline as a ``PulseSchedule``
4. Call ``schedule.time_correct_tally(...)`` to get the dose rate(s)

### Workflow

<!-- doctest: skip -->
```python
import yamc

_chain = "transmutation-endf-b8.1-sfr.arrow"
yamc.transmutation_decay_data = _chain
yamc.transmutation_reactions = _chain
yamc.transmutation_fission_yields = _chain

# 1. Bin a photon tally by parent radionuclide.
#    parent_nuclides=... makes the tally score one bin per radionuclide, so
#    each nuclide's decay photons can be scaled independently in
#    post-processing. The radionuclides come from the model's materials, so
#    no schedule is needed yet. Any spatial filter works: a mesh here, but
#    cells and materials work too.
mesh = yamc.RegularRectangularMesh.from_domain(geometry, shape=1_000_000)
radionuclides = yamc.Model(geometry=geometry, source=source).radionuclides()
dose_tally = yamc.Tally(
    scores=["flux"],
    particle="photon",
    parent_nuclides=radionuclides,
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
flat bin layout per row is score then mesh (the general order is in
[Spatial binning](tallies.md#spatial-binning)), so reshape with
``np.array(dose.mean).reshape(mesh.shape)``.

### Nuclear-data uncertainty on the dose

``dose.std_dev`` is the transport's statistical error alone. Pass
``data_uncertainty`` and the time correction's nuclear-data uncertainty comes
back beside it:

<!-- doctest: skip -->
```python
dose = schedule.time_correct_tally(
    results[dose_tally],
    steps=-1,
    data_uncertainty=yamc.DataUncertainty(seed=1),
)
dose.data_std_dev            # from the nuclear data, same shape as dose.mean
dose.total_std_dev           # std_dev and data_std_dev in quadrature
dose.data_uncertainty_info   # what was sampled, what states no sigma, the
                             # replica count, the convergence target and
                             # whether every factor's sigma met it
```

Two sources act here. ``"half_life"`` acts through the time correction: a
time-correction factor is an activity over the schedule, and the tally's photon
yield is per decay, so the half-lives enter there and nowhere else. Each
replica draws every half-life feeding an emitter once and uses it for every
campaign. ``"decay_photon_lines"`` acts through each photon spectrum's
normalisation (FD for lines, FC for a continuum), which scales an emitter's whole
spectrum and so its tally, and is applied after the fact. The draws are those a
transmutation with the same seed makes. ``dose.data_std_dev_correlated`` and
``dose.total_std_dev_correlated`` give the upper end of the range left by the
unstated correlation between an emitter's gamma and x-ray normalisations. The
line intensities and energies change the spectrum's shape, which would need
line-resolved tallies, and are held; ENDF/B-VIII.1 folds the normalisation's
sigma into the line intensities, so under it almost nothing is drawn. Decay
branching ratios also shape a time correction but are held at nominal here, and
``data_uncertainty_info`` lists everything held under ``not_perturbed``. The two sigmas are
independent (one comes from transport sampling, the other from the evaluated
half-lives), so they add in quadrature. Without ``data_uncertainty``
the three attributes are ``None``.

### Mixed-spectrum (DD + DT) campaigns

A schedule can mix sources, for example a DT pulse, cooling, then a DD pulse. Each
distinct spectrum activates the materials differently, so run transport once
per source and pass ``time_correct_tally`` a dict mapping each source in
``schedule.sources`` to its run's tally result; the per-source corrections are
summed exactly.

<!-- doctest: skip -->
```python
dt = yamc.NeutronSource(energy=yamc.sources.fusion_neutron_spectrum(ion_temp=20000.0, reactants="DT"))
dd = yamc.NeutronSource(energy=yamc.sources.fusion_neutron_spectrum(ion_temp=20000.0, reactants="DD"))

schedule = yamc.PulseSchedule([
    yamc.Pulse(source=dt, rate=1e18, duration=(1, "a")),
    yamc.Cooldown(duration=(1, "d")),
    yamc.Pulse(source=dd, rate=5e17, duration=(30, "d")),
    yamc.Cooldown(duration=(7, "d")),
])

# One transport run per source, the same model otherwise.
model_dt = yamc.Model(geometry=geometry, tallies=[dose_tally], source=dt,
                      transport_secondary_photons=True, use_decay_photons=True)
model_dd = yamc.Model(geometry=geometry, tallies=[dose_tally], source=dd,
                      transport_secondary_photons=True, use_decay_photons=True)
results_dt = model_dt.simulate_transport(total_particles=1_000_000)
results_dd = model_dd.simulate_transport(total_particles=1_000_000)
dose = schedule.time_correct_tally({dt: results_dt[dose_tally],
                                    dd: results_dd[dose_tally]})
```

### A machine that changes between steps

The workflow above activates and transports in one geometry. A real machine
changes: it is irradiated with DD neutrons in one configuration and with DT in
another, components are added between campaigns, and during shutdown, maintenance
and decommissioning they are removed or moved. ``yamc.simulate_shutdown_dose``
takes a schedule whose steps name the ``Model`` the machine is in, and does the
transport itself.

- A ``Pulse(model=...)`` irradiates that model with its own neutron source (give
  ``model`` instead of ``source``).
- A ``Cooldown(model=...)`` is the machine the decay photons are transported
  through while it cools. It needs no source.
- A step without a ``model`` stays in the previous step's.
- Cells are the same component across models when they have the same **name**,
  so name every cell that activates.
- A cell missing from a model has been removed: it emits nothing there, but
  its activation keeps decaying, and counts again if it comes back.
- A cell somewhere else is described by the step's ``moves``,
  ``{cell or name: (dx, dy, dz)}`` in cm or a 4x4 rigid transformation matrix,
  relative to where the cell was built. The model must hold the cell at its moved
  position; each moved decay photon is checked to land inside it.

<!-- doctest: skip -->
```python
dd_phase = yamc.Model(geometry=reactor_without_plug, source=dd_source)
dt_phase = yamc.Model(geometry=reactor_with_plug, source=dt_source)
maintenance = yamc.Model(geometry=reactor_with_plug_pulled_out)  # no source needed
decommissioned = yamc.Model(geometry=reactor_without_plug)

schedule = yamc.PulseSchedule([
    yamc.Pulse(model=dd_phase, rate=1e16, duration=(30, "d")),
    yamc.Pulse(model=dt_phase, rate=1e17, duration=(30, "d")),
    yamc.Cooldown(duration=(1, "d")),                     # still the DT layout
    yamc.Cooldown(model=maintenance, duration=(7, "d"),
                  moves={"port plug": (150, 0, 0)}),       # plug pulled 150 cm out
    yamc.Cooldown(model=decommissioned, duration=(30, "d")),
])

# A photon tally, without parent_nuclides: it is binned by every radionuclide
# the irradiations make. Each tally needs a name; results are keyed by it.
worker = yamc.Tally(scores=["flux"], name="worker", mesh=mesh, particle="photon")
dose = yamc.simulate_shutdown_dose(schedule, tallies=[worker], total_particles=1_000_000)
dose["worker"].mean        # one row per schedule step
dose["worker"].std_dev     # the whole calculation's standard error
dose["worker"].by_nuclide
dose["worker"].times
```

For every model a pulse irradiates in, the neutrons are transported once with
decay photons on, and each decay photon is recorded where it is born instead of
being tracked. Then, for each configuration the machine is observed in, each
neutron history's decay photons, from every irradiated model, are moved with
their cells (or dropped, for removed cells) and transported through the
observing model together, as one photon history. That history's dose at the end
of each step is formed from the time correction factors of each irradiating
configuration's own pulses. A schedule with two irradiating and four observing
configurations takes two neutron runs and four photon runs, however many steps
it has, and every photon is transported once.

``std_dev`` is the standard error of the whole calculation. It is the spread of
the per-history doses, so it covers the neutron transport as well as the photon
transport, and the correlation one neutron history carries between
radionuclides and between irradiations. For a machine that never changes it
equals the error of the single-model run above. There is nothing to tune: the
photon runs have one history per neutron history. Runs on the CPU, on one
process.

### Reference

- ``model.radionuclides()``: activation radionuclides reachable from a Model's materials; pass to ``Tally(parent_nuclides=...)``.
- ``yamc.Pulse(rate, duration, source=None, model=None, moves=None)``: an irradiation campaign. ``duration`` is seconds or a ``(value, unit)`` tuple (units: s, min, h, d, a). ``source`` is the NeutronSource for transport and SDR. ``Model.simulate_transmutation`` may omit it and transports the model's source; ``Material.transmute`` needs one whose energy is a ``Histogram``, and there ``rate`` is the total flux in n/cm^2/s.
- ``yamc.Cooldown(duration, model=None, moves=None)``: a cooling step (no source). ``model`` and ``moves`` describe the machine during the step, for ``simulate_shutdown_dose``.
- ``yamc.simulate_shutdown_dose(schedule, tallies, total_particles, seed=1, threads=None, method="d1s", mesh=None)``: the dose over a schedule whose steps name models; returns ``{tally name: DoseResult}`` with one row per step. ``method="d1s"`` (the default) is the method described above; ``method="r2s"`` solves cell inventories instead, by whole cell or, with a ``RegularRectangularMesh`` as ``mesh``, by cell part under each voxel.
- ``yamc.PulseSchedule([steps...])``: the timeline; ``.sources`` lists its distinct sources, ``len(schedule)`` its step count.
- ``schedule.time_correct_tally(results, steps=None, data_uncertainty=None)``: returns a ``DoseResult`` with ``.mean`` / ``.std_dev`` / ``.by_nuclide`` / ``.times``, and ``.data_std_dev`` / ``.total_std_dev`` / ``.data_uncertainty_info`` when ``data_uncertainty`` is given.

Parent-nuclide binning is set by the ``parent_nuclides=`` argument on
``Tally`` (step 1). No separate setup call is needed.
