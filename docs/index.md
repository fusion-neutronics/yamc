# YAMC

YAMC (Yet Another Monte Carlo) is an open-source Monte Carlo particle transport code for
neutron and photon simulations, aimed at fusion neutronics workflows.

It supports constructive solid and mesh-based geometries, coupled neutron-photon transport,
transmutation, shutdown dose rate analysis, CAD geometry import, flexible tallies, and
convenient post-processing, all accessible through a Python API.

YAMC's priorities are:

- First-class support for fusion neutronics workflows
- Runs on Linux, macOS, and Windows, on x86 and ARM CPUs and AMD, NVIDIA, and Intel GPUs, and in the browser via WebAssembly
- Straightforward installation with `pip`, no compilation step
- A fully permissively licensed open-source dependency stack
- Batteries included: CAD geometry meshing, transmutation, convenient cross section data downloading, and visualization built in

## How these priorities are met

Four design choices make this possible.

The core is written in Rust. All of the physics, geometry, and data handling live in a
Rust workspace. Rust compiles to self-contained native binaries with no runtime or shared-library
dependencies, so the same code base packages cleanly across Linux, macOS, and Windows on both x86
and ARM. The dependency-free core crates also compile to WebAssembly, so the same source can
run in a browser.

The Python package is a thin set of [PyO3](https://pyo3.rs) bindings over that core. The calls are
wrappers and all computation happens in Rust at native speed. The bindings target
Python's stable ABI (abi3, 3.10+), so a single wheel per platform installs with `pip` without
compiling anything.

Apache Arrow is the common data format. Nuclear data, tally results, and meshes are all stored and exchanged as Arrow. Arrow is columnar
and language-agnostic, so data files are portable and load quickly, and tally results convert to
pandas, polars, or xarray with little or no copying.

GPU execution runs through cubecl/Vulkan instead of a proprietary single-vendor toolkit. One code
path therefore targets AMD, NVIDIA, and Intel hardware, and the whole stack stays open source.

## Explore

- [Installation](installation.md): install with `pip`, no compilation step
- [Quick start](quick_start.md): a minimal simulation end to end
- [Usage guide](usage_guide.md): materials, geometry, sources, tallies, transmutation, CAD
- [Nuclear data libraries](nuclear_data_libraries.md): the six you can choose between, and how much is in each
- [Python API](api.md): the type-annotated reference, generated from the package's stubs

To ask questions or get support, try the
[GitHub Discussions](https://github.com/fusion-neutronics/yamc/discussions) page.
