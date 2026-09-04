# Quick start guide

A minimal simulation end to end: configure data, build a material, geometry and source, add a
tally, run, and read the result. The example computes tritium production (the basis of the
**tritium breeding ratio, TBR**) in a lithium sphere.

## Import YAMC

```python
import yamc
```

## Configure the nuclear data

A single library keyword downloads the required data automatically:

```python
yamc.cross_section_data = "endf-b8.1"
```

There are six to choose between, and they differ by more than a name:
[Nuclear data libraries](nuclear_data_libraries.md).

## Create a material

Lithium breeds tritium, so use it as the breeder material (`"Li"` expands to natural Li6 + Li7):

```python
mat_lithium = yamc.Material(
    composition={"Li": 1.0},
    density=0.534,
    temperature=294,
)
```

## Build the geometry and source

```python
surface_sphere = yamc.Sphere(radius=100, boundary="vacuum")
cell = yamc.Cell(region=surface_sphere.below, material=mat_lithium)
geometry = yamc.Geometry(cells=[cell])

source = yamc.NeutronSource(energy=14.06e6)  # a float is a monoenergetic source (14.06 MeV, D-T)
```

## Add a tally

Score H3 (tritium) production in the lithium cell. Per source neutron, this is the TBR:

```python
tally = yamc.Tally(scores=["H3-production"], cells=cell)
```

## Run the simulation

Pass the tally to the model and run the transport:

```python
model = yamc.Model(
    geometry=geometry,
    source=source,
    tallies=[tally],
)
results = model.simulate_transport(total_particles=100_000)
```

## Read the result

`results[tally]` is a labelled result. Its `.mean` is the H3-production rate per source neutron
(the TBR for this geometry) and `.standard_deviation` the statistical uncertainty. Each is a
list with one entry per score:

```python
result = results[tally]
print(f"TBR (H3 production per source neutron): {result.mean} +/- {result.standard_deviation}")
```

For a pure-lithium sphere this prints a TBR of roughly:

```text
TBR (H3 production per source neutron): 1.223 +/- 0.014
```
