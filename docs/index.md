# YAMC

**YAMC** (Yet Another Monte Carlo) is an open-source Monte Carlo particle transport code for
**neutron and photon simulations**, with a focus on **fusion neutronics workflows**.

It supports constructive solid and mesh-based geometries, coupled neutron-photon transport,
transmutation, shutdown dose rate analysis, CAD geometry import, flexible tallies, and
convenient post-processing, all accessible through a Python API.

YAMC's priorities are:

- First-class support for fusion neutronics workflows
- Runs everywhere: **Linux, macOS, and Windows**, on x86 and **ARM** CPUs and AMD, NVIDIA, and Intel **GPUs**, and in **browser** via WebAssembly
- Straightforward installation with `pip`, no compilation step
- A fully permissively open-source dependency stack

- Batteries included: CAD geometry meshing, transmutation, convenient cross section data downloading, and visualization built in

## How these priorities are met

A few deliberate design choices let YAMC meet all of the above at once:

**A Rust core, compiled everywhere.** All of the physics, geometry, and data handling live in a
Rust workspace. Rust compiles to self-contained native binaries with no runtime or shared-library
dependencies, so the same code base packages cleanly across Linux, macOS, and Windows on both x86
and ARM. The dependency-free core crates also compile to WebAssembly, which is what makes
browser-based deployment possible from the very same source.

**A thin Python API over that core.** The Python package is a set of [PyO3](https://pyo3.rs)
bindings: the calls are thin wrappers and all computation happens in Rust, so you get native
performance through a familiar Python neutronics workflow. The bindings target Python's stable ABI
(abi3, 3.10+), so a single wheel per platform installs with `pip` and no compilation step, keeping
installation straightforward.

**Apache Arrow as the common data format.** Nuclear data, tally results, and meshes are all stored
and exchanged as Arrow. Because Arrow is columnar and language-agnostic, data files are portable
and load quickly, and tally results convert to pandas, polars, or xarray with little or no copying
-- so you can stay in whichever analysis tools you already use.

**An open dependency stack.** GPU execution runs through cubecl/Vulkan rather than a proprietary,
single-vendor toolkit, which is what allows one code path to target AMD, NVIDIA, and Intel hardware
while keeping the whole stack open source.

## Explore

- [Installation](installation.md): install with `pip`, no compilation step
- [Quick start](quick_start.md): a minimal simulation end to end
- [Usage guide](usage_guide.md): materials, geometry, sources, tallies, transmutation, CAD
- [Nuclear data libraries](nuclear_data_libraries.md): the six you can choose between, and how much is in each
- [Python API](api.md): the type-annotated reference, generated from the package's stubs

To ask questions or get support, try the
[GitHub Discussions](https://github.com/fusion-neutronics/yamc/discussions) page.
