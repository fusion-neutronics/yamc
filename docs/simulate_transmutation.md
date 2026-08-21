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

### Per-edge reaction rates

The solve computes the rate of every individual production edge to build its burnup matrix, whose entries are sums over the edges into each product. ``get_reaction_rates`` hands those edges back, so a route can be weighted rather than merely listed:

<!-- doctest: skip -->
```python
rates = results.get_reaction_rates(material_id=1, step=0)   # step indexes timesteps
rates["Fe58"]["(n,2n)"]        # [('Fe57', 1.205e-09)]      # rate in 1/s
```

``parent -> reaction kind -> [(target, rate [1/s])]``. A ``target`` of ``None`` is a channel naming no single product, which in practice means fission, whose products come from the yields rather than from an edge. A decay-only step drove no reaction and reports an empty map.

Note the index: ``step`` is the schedule step, the same index as ``timesteps`` and ``source_rates``, which is one less than the ``step`` the composition getters take (where 0 is the initial composition).

Both methods report the same quantity in the same shape. ``method="coupled"`` takes the rates from that step's own transport tallies and ``method="independent"`` from the single transport's multigroup fold scaled by the step's source rate. Either way an edge rate is that rate times the branching of the chain the step was solved with, isomeric overlay included, so a consumer does not have to know which method produced the result.

The share of a product arriving down each route is the headline column of a decay heat pathway analysis, and it is what makes a C/E actionable: knowing a foil reads 2.5x the measurement says something is wrong, and knowing that 73% of the dominant product arrives down one channel says which evaluated quantity to go and look at.

<!-- doctest: skip -->
```python
from collections import defaultdict

into = defaultdict(list)
for parent, kinds in rates.items():
    for kind, edges in kinds.items():
        for target, rate in edges:
            if target is not None:
                into[target].append((rate, parent, kind))

routes = sorted(into["Fe57"], reverse=True)
total = sum(rate for rate, _, _ in routes)
for rate, parent, kind in routes:
    print(f"{rate / total:6.1%}  {parent}{kind}Fe57   {rate:.4g} /s")
```

For a 14.06 MeV point source in a 10 cm iron sphere, one hour at 1e18 n/s:

```text
 73.3%  Fe58(n,2n)Fe57       1.205e-09 /s
 18.3%  Co58(n,np)Fe57       3.016e-10 /s
  6.1%  Ni60(n,a)Fe57        1.009e-10 /s
  1.3%  Co58(n,d)Fe57        2.176e-11 /s
  0.6%  Fe56(n,gamma)Fe57    9.747e-12 /s
  0.3%  Ni61(n,na)Fe57       4.809e-12 /s
  0.0%  Co59(n,t)Fe57        6.359e-13 /s
  0.0%  Ni59(n,3He)Fe57      8.149e-15 /s
```

A rate is per atom of the parent, so weighting a route by inventory means multiplying by the parent's density at that step (``get_nuclide_density``). Over one short step the two together reproduce the step's own burnup: the Fe57 produced is the summed edge rate times the parent densities times the step length.

### Material.transmute()

For standalone transmutation with a pre-computed flux spectrum (no transport), see [Material transmutation](materials.md#transmutation) in the Materials section.
