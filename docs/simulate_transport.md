
## Running the Simulation

### model.simulate_transport()

<!-- doctest: skip -->
```python
results = model.simulate_transport(total_particles=10000)
```

``simulate_transport()`` returns a ``SimulationResults`` object containing a ``TallyResult`` for every tally attached to the model, plus run metadata (``n_batches``, ``particles_per_chunk``, ``elapsed``).

Every run needs at least one **stop condition**: a particle budget (``total_particles``), a wall-time budget (``max_runtime``), or convergence targets on the model. They are OR-combined, so the run ends at the first one satisfied. ``total_particles`` defaults to ``None`` (no particle cap); a call with no stop condition at all raises ``ValueError``.

#### Accessing tally results

A ``TallyResult`` can be looked up four ways: by the ``Tally`` object, by name, by id, or by iteration. Names and ids are optional -- set them on the ``Tally`` constructor to use those flavors.

<!-- doctest: skip -->
```python
flux_tally = yamc.Tally(scores=["flux"], cells=cell, name="flux", id=1)
heating_tally = yamc.Tally(scores=["heating"], cells=cell)

model = yamc.Model(
    geometry=geometry,
    source=source,
    tallies=[flux_tally, heating_tally],
)
results = model.simulate_transport(total_particles=10000)

# By Tally object -- always works
flux_result = results[flux_tally]

# By name -- requires name= on the Tally
flux_result = results["flux"]

# By id -- requires id= on the Tally
flux_result = results[1]

# Iterate over (tally, result) pairs
for tally, result in results.items():
    print(tally.name, result.mean)

# Membership test
assert flux_tally in results
assert "flux" in results
```

A ``TallyResult`` exposes ``mean``, ``standard_deviation``, ``relative_error``, ``variance``, ``shape``, ``dim_labels``, and ``to_numpy()``, plus a set of statistical-reliability metrics and a PASS/FAIL verdict that print whenever you display a result (see [Statistical reliability](results.md#statistical-reliability) and [Results](results.md) below).

### Progress output (``verbose``)

The ``verbose`` argument on ``Model`` is a **list of independent flags** selecting which categories of progress are printed during transport. Progress is reported in **source particles** (the run is internally chunked, but you never deal with batches). The flags:

- ``"progress"`` -- periodic ``Progress: N/total particles (P%)`` lines.
- ``"eta"`` -- progress lines that also show elapsed wall-clock and ETA.
- ``"tally"`` -- each tally's running mean ± std (rel. err) as it converges.
- ``"summary"`` -- the end-of-run completion line and per-tally final statistics.

Mix and match freely. The default (pass nothing, or ``None``) is ``["eta", "tally", "summary"]``. An **empty list ``[]`` is fully silent** -- not even the end-of-run summary is printed.

<!-- doctest: skip -->
```python
# Default output
model = yamc.Model(geometry=geometry, source=source, tallies=[tally])

# Fully silent -- useful for scripts that post-process results themselves
model = yamc.Model(geometry=geometry, source=source, tallies=[tally], verbose=[])

# Just the final summary, no per-particle progress
model = yamc.Model(geometry=geometry, source=source, tallies=[tally], verbose=["summary"])

# Read or change after construction
print(model.verbose)                       # ['eta', 'tally', 'summary']
model.verbose = ["progress", "summary"]
```

### Threading

By default all available CPU cores are used. Specify a thread count explicitly:

<!-- doctest: skip -->
```python
model.simulate_transport(total_particles=10000, threads=8)
```

### GPU transport (``compute="gpu"``)

The transport loop has a GPU implementation that runs on Vulkan via cubecl. Selected with ``compute="gpu"`` on ``simulate_transport``:

<!-- doctest: skip -->
```python
# Check that a usable GPU is present at module level.
if yamc.parallel.gpu_available():
    results = model.simulate_transport(total_particles=1_000_000, compute="gpu")
else:
    results = model.simulate_transport(total_particles=1_000_000)
```

On a machine with more than one GPU -- e.g. an integrated GPU plus a dedicated one -- ``compute="gpu"`` auto-selects an adapter (preferring a discrete GPU via the high-performance power preference). To pin a specific adapter, list the available ones with ``yamc.parallel.list_gpu_adapters()`` and pass its name directly as ``compute=<name>``:

<!-- doctest: skip -->
```python
adapters = yamc.parallel.list_gpu_adapters()
print(adapters)
# e.g. ['NVIDIA GeForce RTX 4070', 'Intel(R) UHD Graphics']

# Pin a specific GPU by name (e.g. the dedicated one):
results = model.simulate_transport(total_particles=1_000_000, compute=adapters[0])
```

``yamc.parallel.gpu_available()`` returns ``True`` only when a GPU with **f64 compute** support is reachable. The kernel uses double-precision throughout (cross-section interpolation, sampling, and the variance accumulator), so single-precision-only consumer GPUs are rejected.

The GPU path supports a defined slice of the option space, **verified to match the CPU path** on real discrete hardware (NVIDIA RTX A500): `cpu_gpu_equivalence_two_materials` is bit-exact and the `gpu_*.rs` GPU-vs-CPU suite passes. As of today it supports:

- **Neutron transport** and **standalone primary-photon transport** (a pure `PhotonSource`). Coupled neutron→photon secondary production (`transport_secondary_photons=True`) and D1S decay photons (`use_decay_photons=True`) are **not** supported - the GPU has no particle bank - and are rejected.
- **Tallies** with the **track-length** estimator: flux, heating / heating-local, damage-energy, reaction rates (total / elastic / absorption / (n,γ) / inelastic / fission), gas-production, and photon cross-section scores. The `collision` estimator is rejected.
- **Tally filters**: `cells=`, `materials=`, `mesh=`, `energy_bins=`, `particle=`, `nuclides=`, `parent_nuclides=`, and `energy_function=` / `dose_coefficients=`. The energy function weights rather than bins, so it still needs a spatial binner alongside it; its cubic spline is solved on the host and only evaluated in the kernel, so GPU and CPU interpolate identically.
- Full neutron physics: elastic (tabulated angular), inelastic (Kalbach-Mann, n-body, evaporation), free-gas thermal, URR, and the photon interactions on the photon path.
- CSG geometry only (mesh is rejected); surfaces sphere, plane, cylinder, Z-torus.
- Analog transport only - `variance_reduction` is rejected (not silently ignored).

If any unsupported option is used, ``compute="gpu"`` raises a ``ValueError`` at simulate-transport start with a specific message. If no GPU with f64 compute is available, it raises ``RuntimeError``. The ``threads`` and ``capture_tracks`` kwargs are not used with ``compute="gpu"`` (the kernel parallelises across all particles in one launch, and tracking is CPU-only). Repeated `compute="gpu"` calls in one process are safe - each resets the shared GPU client.

For the full, always-current per-option breakdown of which CPU options the GPU path supports (and which it rejects), see the CPU-vs-GPU option matrix harness `crates/yamc/tests/gpu_cpu_comparison_matrix.rs`, which is the source of truth.

### MPI (distributed memory)

For multi-node parallelism, build YAMC with MPI support (see [installation](installation.md)).

<!-- doctest: skip -->
```python
# Check MPI support
if yamc.parallel.has_mpi():
    print(f"Rank {yamc.parallel.mpi_rank()} of {yamc.parallel.mpi_size()}")
```

Run with ``mpirun``:

```bash
mpirun -np 4 python simulation.py
```

Hybrid MPI + threading:

<!-- doctest: skip -->
```python
model.simulate_transport(total_particles=10000, threads=8)  # 8 threads per MPI rank
```

```bash
# 4 ranks x 8 threads = 32 cores
mpirun -np 4 python simulation.py
```

Tallies are automatically reduced to rank 0 after each batch. Always check ``yamc.parallel.mpi_rank() == 0`` before reading results.

Results are deterministic regardless of rank count -- each particle gets a unique seed based on the base seed, batch number, and global particle index.

Slurm example:

```bash
#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=32
srun python simulation.py  # In script: model.simulate_transport(total_particles=1_000_000, threads=32)
```

### Particle tracking

Track particle histories for debugging or visualization. Tracks are off by default (zero overhead) and are only captured when ``capture_tracks=`` is set on ``simulate_transport()``. Histories are addressed by their global 0-based index, so the captured set is identical regardless of thread count:

- ``capture_tracks=N`` -- the first ``N`` histories (recommended for debugging)
- ``capture_tracks=range(20, 30)`` -- histories 20..29 (half-open); strides work too, e.g. ``range(0, 1000, 10)`` for every 10th of the first 1000
- ``capture_tracks=[20, 30, 45]`` -- exactly those history indices
- ``capture_tracks='all'`` -- every history (warns; high memory on large runs)
- ``capture_tracks=None`` / omitted -- no tracking

When tracking is enabled the captured ``Tracks`` are attached to the returned ``SimulationResults`` as ``results.tracks``; otherwise ``results.tracks`` is ``None``.

<!-- doctest: skip -->
```python
results = model.simulate_transport(total_particles=10000, capture_tracks=100)
tracks = results.tracks  # None if capture_tracks= was not set

print(len(tracks))            # number of particle histories captured
print(tracks.total_events())  # total events across all histories

# Iterate over individual particle tracks
for track in tracks.tracks[:3]:  # first 3 particles
    print(f"Particle {track.particle_id} (generation {track.generation}):")
    for event in track.events[:5]:  # first 5 events
        print(f"  {event.event_type:12s}  E={event.energy_in:.4e} -> {event.energy_out:.4e} eV"
              f"  pos=({event.x:.2f}, {event.y:.2f}, {event.z:.2f})"
              f"  cell={event.cell_id}  nuclide={event.nuclide}  MT={event.reaction_mt}")

# Convert to records for pandas
import pandas as pd
df = pd.DataFrame(tracks.to_dataframe_records())
print(df[["event_type", "energy_in", "energy_out", "x", "y", "z", "nuclide"]].head(10))
```

Each ``TrackEvent`` includes: ``event_type``, ``position``, ``direction``, ``energy_in``, ``energy_out``, ``weight``, ``cell_id``, ``reaction_mt``, ``nuclide``, and ``birth_reaction``.

### Lost particles and geometry debugging

If a particle crosses into a region not covered by any cell, it is "lost". YAMC prints diagnostics and kills the particle.

<!-- doctest: skip -->
```python
model = yamc.Model(
    geometry=geometry,
    source=source,
    max_lost_particles=10,  # abort after this many (default 10)
)
model.simulate_transport(total_particles=1000)

# Inspect lost particles
for lp in model.lost_particles:
    print(f"Position: {lp.position}, Surface: {lp.surface_id}")
```


