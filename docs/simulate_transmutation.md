# Transmutation 

### Model.simulate_transmutation()

Coupled transport-transmutation runs transport at each timestep with updated material compositions. Materials must have ``transmutable=True`` and a ``volume`` set.

<!-- doctest: skip -->
```python
_chain = "transmutation-endf-b8.1-sfr.arrow"
yamc.transmutation_decay_data = _chain
yamc.transmutation_reactions = _chain
yamc.transmutation_fission_yields = _chain

model = yamc.Model(geometry=geometry, source=source)

schedule = yamc.PulseSchedule([
    yamc.Pulse(rate=1e14, duration=(1, "d"), source=source),   # irradiation, rate in n/s
    yamc.Pulse(rate=1e14, duration=(1, "d"), source=source),   # irradiation
    yamc.Cooldown(duration=(1, "d")),                          # cooling
])
results = model.simulate_transmutation(
    method="coupled",                    # or "independent"
    schedule=schedule,
    total_particles=5000,
)
```

The irradiation/cooling timeline is one ``PulseSchedule`` of ``Pulse`` / ``Cooldown`` steps, the same type used for the shutdown-dose-rate post-processing below. Each ``Pulse`` carries a duration (seconds, or a ``(value, unit)`` tuple) and a source ``rate`` in n/s; each ``Cooldown`` is decay-only. Transport uses the source configured on the ``Model``, so a schedule with more than one distinct ``Pulse`` source is rejected.

Methods:
- ``"coupled"`` -- transport runs every step with updated compositions
- ``"independent"`` -- transport runs once; reaction rates scaled by source rate per step (faster for low-burnup / fusion)

### TransmutationResults

<!-- doctest: skip -->
```python
# Evolution of a nuclide over all timesteps
evolution = results.get_nuclide_evolution(material_id=1, nuclide="Co60")
print(results.times)  # cumulative times including t=0

# Density at a specific step
density = results.get_nuclide_density(material_id=1, nuclide="Co60", step=2)

# Full material at a step
mat = results.get_material(material_id=1, step=2)
print(mat.nuclides)

# All nuclides at a step
nuclides = results.get_material_nuclides(material_id=1, step=2)
```

### Material.transmute()

For standalone transmutation with a pre-computed flux spectrum (no transport), see [Material transmutation](materials.md#transmutation) in the Materials section.
