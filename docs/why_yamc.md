# Why YAMC

YAMC is a Monte Carlo neutron and photon transport code for fusion. These are
the things it does well.

## Easy

**Try it in seconds, with nothing installed.**
[YAMC online](https://fusion-neutronics.github.io/yamc-online/) runs the same
transport engine in your browser, on your own machine, with no account and no
upload.

**`pip install yamc` on Linux, macOS or Windows.** The wheel carries the
compiled solver, so there is no compiler or toolchain to set up. See
[Installation](installation.md).

**Nuclear data downloads itself.** Pick a library and only the parts your
model needs are fetched and cached, down to just the temperatures your
materials are at. To work offline, download the full data files and point
YAMC at the directory. See [Nuclear data](nuclear_data.md).

**Permissively licensed throughout.** YAMC is MIT or Apache-2.0, and nothing
in its dependency stack is copyleft, including the GPU path.

**Share a model as a web page.** `model.to_html()` writes one file that runs
the simulation in any browser, for a colleague without Python or a demo in a
course. See
[Export to a runnable HTML page](models.md#export-to-a-runnable-html-page-to_html).

## Accurate

**Continuous-energy neutron and photon transport.** Cross sections are used on
the evaluation's own energy grid, and the photon physics, including atomic
relaxation, is read straight from the ENDF sublibraries with no NJOY step. See
[Method](method.md).

**Every tally says whether to trust it.** Statistics are kept per history, with
no batches to choose, and each result carries automatic reliability checks on
its error, variance of variance and convergence. See
[Statistical reliability](results.md#statistical-reliability).

**Uncertainty from the nuclear data.** A transmutation can return every
nuclide density, activity, decay heat and dose with an uncertainty propagated
from the evaluation's own covariance data. See
[Nuclear-data uncertainty](materials.md#nuclear-data-uncertainty).

## Fast

**GPU transport on any vendor's hardware.** `compute="gpu"` runs on AMD,
NVIDIA or Intel GPUs in double precision. See
[GPU transport](simulate_transport.md#gpu-transport-computegpu).

**Variance reduction built in.** Weight windows are generated automatically
from the model, and delta tracking and survival biasing are one setting each.
See [Weight windows](weight_windows.md).

**Stop when the answer is good enough.** Set a target relative error on a
tally and the run ends as soon as it is met. See
[Stopping on a precision target](results.md#stopping-on-a-precision-target).

**Shutdown dose from one transport run.** The dose at any cooling time, for
any irradiation schedule, comes from a single run. See
[Shutdown dose rate](simulate_shutdown_dose.md).

## Flexible

**CSG, CAD and meshes in one model.** YAMC includes its own CAD mesher and
mesh transport, so a CadQuery assembly or STEP file goes straight into a
model, alongside or inside CSG cells. See [Geometry](geometry.md).

**Built for fusion.** A tokamak plasma source from plasma parameters, and
tallies for heating, DPA, tritium production and dose. See
[Sources](sources.md) and [Tallies](tallies.md).

**Transmutation with or without transport.** Deplete materials coupled to a
transport run, or from a spectrum you already have, and get activity, decay
heat, decay photons and waste clearance. See
[Simulate transmutation](simulate_transmutation.md).

**Mix nuclear data libraries.** Cross sections, reactions, decay data and
fission yields can each come from a different library, or from your own files.
See [Nuclear data settings](nuclear_data.md#nuclear-data-settings).

**Results where your analysis already is.** Tallies convert to numpy, pandas,
polars or xarray, and export to VTK-HDF for ParaView. See [Results](results.md).
