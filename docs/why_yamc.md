# Why YAMC

Four things, in order of how often they decide the question: it has to be
easy to install and run, the numbers have to be right and come with an honest
error bar, it has to get to that error bar quickly, and it has to bend to the
fusion problem you actually have rather than the one a tool author imagined.

## Easy

**`pip install yamc`, on Linux, macOS or Windows.** The wheel carries the
compiled Rust solver, threading included. There is no compiler to invoke, no
toolchain to build first and no system package to hunt down, and one wheel per
platform covers every Python from 3.10 up. See
[Installation](installation.md).

**Nuclear data downloads on demand.** Set `yamc.cross_section_data =
"endf-b8.1"`, or any of the other five library keywords, and the nuclides your
model contains are fetched and cached on first use, not a whole library up
front. Only the sections a run uses come down: an activation run on Fe56 skips
the 15 MB of secondary distributions that transport needs. See [Nuclear data libraries](nuclear_data_libraries.md).

**The license stack is fully permissive.** YAMC is dual-licensed under MIT or
Apache-2.0, at your option, and nothing in the dependency chain forces a
copyleft license onto whatever you build with it. The GPU path runs on
cubecl/Vulkan rather than a proprietary vendor toolkit, so that holds for GPU
runs too. Using it in a closed pipeline or a commercial tool is not a licensing
question, just an engineering one.

**A citable material library ships in the wheel.** The PNNL Compendium
(PNNL-15870 Rev. 2), 410 named materials with their compositions and densities,
is bundled rather than a separate download, so a lookup needs no network and
cannot fail. A standard shielding or structural material is
`yamc.materials.pnnl["Steel, Stainless 304"]` instead of a retyping exercise.
See [Material collections](materials.md#material-collections).

**A model can be handed to someone who has nothing installed.**
`model.to_html()` writes one self-contained HTML page that runs the simulation
in the browser, on a WebAssembly build of the same transport engine. It is a
working example for a colleague without Python, a runnable demo for a paper or
a course, and with `embed_cross_sections=True` it works offline. See
[Export to a runnable HTML page](models.md#export-to-a-runnable-html-page-to_html).

**A model is a file you can diff.** `model.save()` writes the whole
specification as JSON, so a model can be checked into git and reviewed, built
in a notebook and run on a cluster, or published alongside a paper. Results go
to one Arrow file with `results.to_arrow()`, statistics and provenance
included. See [Save and load](models.md#save-and-load-json).

## Accurate

**Continuous-energy neutron and photon transport.** Cross sections are
interpolated on the evaluation's own grid, with unresolved-resonance
probability tables, free-gas scattering at low energy, and the evaluated
secondary-neutron laws. Summed reactions are rebuilt from their partials using
the ENDF sum rules, so nothing is counted twice. See [Method](method.md).

**Temperatures between library points are blended, not sampled.** A material
between two library temperatures gets a new data set built once at load,
linear in temperature. Picking one of the two bracketing temperatures at random
for each lookup, the cheaper alternative, biases deep-penetration results,
because attenuation through a thick shield is an average of exponentials and
not the exponential of an average. See [Temperature](method.md#temperature).

**Photon physics read straight from the evaluation.** Compton scattering with
Doppler broadening, the full cascade of radiative and non-radiative transitions
after a photoelectric absorption, pair production and thick-target
bremsstrahlung, all from the ENDF photo-atomic and atomic relaxation
sublibraries with no NJOY step. Photon data defaults to ENDF/B-VIII.1 whatever
the neutron library is, so fluorescence still works when the neutron data is
TENDL or FENDL. See [Photon transport](method.md#photon-transport).

**Statistics per history, with no batches to choose.** Every tally bin keeps a
per-history mean and variance, which is correct when secondaries are
correlated with their parent, as fission and (n,xn) progeny are. There is no
batch size to tune and no batch-to-batch estimate to distrust. See
[Mean and standard deviation](method.md#mean-and-standard-deviation).

**Every tally says whether to trust it.** Each result carries statistical
checks modelled on MCNP's: relative error and variance of variance below 0.1,
a tail slope above 3, and a mean, error and figure of merit that are not
drifting over the second half of the run. They are computed from a small
per-history aggregate, so they cost next to nothing and are always there.
`result.statistical_checks()` reports them. See
[Statistical reliability](results.md#statistical-reliability).

**A history is the same whatever runs it.** Each history has its own random
number stream, seeded from the base seed and its index in the run. The thread,
chunk and MPI rank do not enter, so changing the core count or the cluster
layout changes nothing but the last bits of the floating-point sums. The CPU
and GPU kernels share the generator and follow the same histories. See
[Random numbers and reproducibility](method.md#random-numbers-and-reproducibility).

**Transmutation rates from continuous energy, not a group fold.** In a
transport-coupled depletion run the reaction rates are scored as
$w\,\ell\,\sigma(E)$ on every track, and isomeric branching ratios are applied
at the collision energy. Nothing is collapsed onto a group structure first, so
resonance structure inside a group is not averaged away before the rates are
formed. See [Transmutation coupling](method.md#transmutation-coupling).

**An uncertainty on the inventory, not just the inventory.** Pass
`data_uncertainty=` to a transmutation and every nuclide density comes back
with a standard deviation, from the evaluation's own cross-section covariance,
the decay data's own sigmas and, when you give one, the spectrum's own
error. Activity, decay heat and contact dose come back with a band
from the same ensemble. See
[Nuclear-data uncertainty](materials.md#nuclear-data-uncertainty).

**It refuses rather than quietly substitutes.** A spectrum above an
evaluation's energy ceiling, or a transmutation that needs data that has been
switched off, raises an error saying why rather than silently returning a
different calculation, and an unsupported estimator or GPU combination is
refused or documented, never swapped for another. Where a library is known to be incomplete, as with
placeholder decay energies, the docs say how much and where. See
[What a library can quietly lack](nuclear_data_libraries.md#what-a-library-can-quietly-lack).

## Fast

**GPUs from any vendor, in double precision.** `compute="gpu"` runs the
transport loop on AMD, NVIDIA or Intel hardware through Vulkan, with
double precision throughout, and it is verified against the CPU path on
discrete hardware. It covers a defined subset of the options and rejects
the rest with a specific message. See
[GPU transport](simulate_transport.md#gpu-transport-computegpu).

**Delta tracking where it pays, surface tracking where it does not.**
`tracking_mode="hybrid"` delta-tracks the cells where it is cheap and
surface-tracks the ones where it is not, which is the safe way to get
delta-tracking speed on a realistic fusion model. `"woodcock"` can be several
times faster than surface tracking on dense, void-free models such as a finely
diced lattice. All three modes give the same answer within statistics.
See [Tracking mode](models.md#tracking-mode-woodcock-delta-tracking).

**Weight windows generated for you.** `model.generate_weight_windows()` builds
neutron and photon windows on a mesh from runs at reduced density
(DeGVR), choosing the density reduction from the model's own optical depth, so
a deep-penetration shielding problem does not need a hand-tuned importance map
or a separate deterministic solve. Survival biasing is one more line, and on
absorption-dominated problems it typically halves the CPU time for the same
error bars. See [Weight windows](weight_windows.md) and
[Survival biasing](models.md#survival-biasing-implicit-capture).

**Stop when the answer is good enough.** A `ConvergenceTarget` on a tally's
relative error, standard deviation or variance of variance ends the run as
soon as it is met, alongside a particle count or a wall-clock limit, on CPU,
GPU or MPI. See
[Stopping on a precision target](results.md#stopping-on-a-precision-target).

**One transport run for any cooling schedule.** Shutdown dose uses a D1S
method that computes the time correction factors for every reachable
radionuclide after the run, so the irradiation and cooling schedule can change,
and the dose can be had at any number of cooling times, without re-running
transport. See [Shutdown dose rate](simulate_shutdown_dose.md).

**A cluster without MPI.** `yamc.combine_results()` pools separate runs into
one result that is statistically identical to a single longer run, so N jobs
with seeds 1 to N, saved with `results.to_arrow()` and combined afterwards, are
a cluster run with no MPI build at all. It refuses runs that cannot be pooled:
the same seed twice, or a different model or library. See
[Combining results](models.md#combining-results-from-separate-runs).

**Nuclear data is a binary columnar format, not text to parse.** Cross
sections, transmutation chains, meshes and tally results are Apache Arrow, so
loading data is reading typed columns into memory, and only the sections a run
needs are read. Converting from ENDF is a one-time step.

## Flexible

**CSG, CAD and meshes in one model.** Build geometry from planes, spheres,
cylinders, cones, tori and general quadrics, or mesh a CadQuery assembly or a
STEP file with the built-in mesher and track on the surface mesh directly, with
no separate DAGMC build. A mesh can also fill a CSG cell. See
[Geometry](geometry.md) and [CAD geometry](geometry.md#cad-geometry).

**A tokamak plasma source from the plasma parameters.**
`yamc.sources.tokamak_source()` takes the major and minor radius, elongation,
triangularity, Shafranov shift, confinement mode and density and temperature
profiles, and builds the D-T and D-D emission with a Ballabio spectrum in each
voxel. A toroidal sector of the machine gets the matching wedge of plasma. See
[Parametric tokamak plasma source](sources.md#parametric-tokamak-plasma-source).

**Tallies for the quantities a fusion design needs.** Flux, heating, DPA, gas
and tritium production and any MT number, on cells, materials, regular,
cylindrical and unstructured meshes, binned by energy, nuclide or particle.
ICRP dose coefficients and arbitrary response functions weight a tally
directly, and `response=` scores a material's reaction rate in a void, as a
virtual detector. See [Tallies](tallies.md).

**Group structures that match OpenMC's.** Twelve named group structures are
built in, and every one except `CCFE-24-PHOTON` has the same boundaries as
OpenMC's entry of that name, so a spectrum tabulated for one code can be
handed to the other without re-binning. See
[The built-in group structures](tallies.md#the-built-in-group-structures).

**Transmutation with or without transport.**
`model.simulate_transmutation()` couples depletion to transport, re-running
transport each step or once for a low-burnup fusion case.
`material.transmute()` skips transport entirely and takes a spectrum you
already have. Either way the result gives activity, decay heat, decay photon
lines and clearance indices against UK, German, US, EU and IAEA limits. See
[Simulate transmutation](simulate_transmutation.md) and
[Clearance](materials.md#clearance).

**Mix nuclear data libraries by section, and by nuclide.** Cross sections and
the four transmutation subsections are separate settings, so a TENDL reactions
network can borrow ENDF/B-VIII.1 decay data, and `cross_section_data` can be
set per nuclide, mixing downloaded keywords with your own converted files. See
[Nuclear data settings](nuclear_data.md#nuclear-data-settings).

**Results go where your analysis already is.** Tally results convert to numpy,
pandas, polars or xarray, and geometry, sources and mesh tallies export to
VTK-HDF for ParaView. Geometry, models and mesh tallies also plot
interactively in a notebook. See [Results](results.md).
