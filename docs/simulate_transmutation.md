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

### Production routes from the chain

``TransmutationChain`` carries the topology a solve steps through: which nuclide becomes which, and how. ``reactions`` holds the neutron-induced edges, ``decays`` the radioactive ones, both keyed by parent and both giving ``(kind, target, branching)``:

<!-- doctest: skip -->
```python
chain = yamc.TransmutationChain("transmutation-endf-b8.1-sfr.arrow")

chain.reactions["W186"]   # [('(n,2n)', 'W185', 1.0), ..., ('(n,a)', 'Hf183', 1.0)]
chain.decays["Hf183"]     # [('beta-', 'Ta183', 1.0)]
chain.decays["W185_m1"]   # [('IT', 'W185', 1.0)]
```

Stable nuclides have no decay modes and are absent from ``decays``, as they are from ``half_lives``. A ``target`` of ``None`` means the channel names no single product: fission on the reaction side, a product outside the chain on the decay side. Decay modes carry the evaluation's own spellings (``"beta-"``, ``"ec/beta+"``, ``"alpha"``, ``"IT"``, ``"sf"``, and emissions written as ``"beta-,n"``), and the branchings for one parent sum to 1.

Because the two share a shape they can go in one edge list, which is what a route walk wants. With reactions alone a walk stops at the first product, one step short of most real routes:

<!-- doctest: skip -->
```python
edges = {}
for parent, es in chain.reactions.items():
    edges.setdefault(parent, []).extend(es)
for parent, es in chain.decays.items():
    edges.setdefault(parent, []).extend(es)

for kind1, mid, _ in edges["W186"]:              # two-step routes from W186 to Ta183
    for kind2, end, _ in edges.get(mid, []):
        if end == "Ta183":
            print(f"W186{kind1}{mid}({kind2}){end}")

# W186(n,3n)W184((n,np))Ta183
# W186(n,4n)W183((n,p))Ta183
# W186(n,a)Hf183(beta-)Ta183      <- the one that needs the decay edge
```

The branching on an edge is a share of its channel, not a rate: it says how a channel splits, not how often the channel fires. Two routes into the same product are not weighted against each other by these numbers alone.

### Material.transmute()

For standalone transmutation with a pre-computed flux spectrum (no transport), see [Material transmutation](materials.md#transmutation) in the Materials section.
