# WebAssembly

yamc has WebAssembly (WASM) bindings, so the library can run in web browsers and other JavaScript environments. Nuclear data processing and visualization then run on the client, with no server-side computation.

## Features

The WebAssembly version supports:
- Material creation and manipulation
- Cross section calculation and visualization
- Predefined materials
- Interactive plotting with Plotly
- Loading an exported model and running the transport loop in the browser

Not every Rust API is exposed to WebAssembly; the bindings cover the
material, plotting, and in-browser simulation paths described below.

## Building WebAssembly

The WebAssembly bindings are compiled from the Rust in
[fusion-neutronics/core](https://github.com/fusion-neutronics/core), so the
commands below run in a checkout of that repository, not of this one:

```bash
git clone https://github.com/fusion-neutronics/core.git
cd core
```

### Prerequisites

Install wasm-pack:

```bash
curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
```

### Building the Package

From the checkout's `crates/yamc` directory:

```bash
cd crates/yamc
# A quick, unoptimised build
wasm-pack build --dev --target web --features wasm
# The optimised build: wasm-pack builds in release mode unless told otherwise
wasm-pack build --target web --features wasm
```

Keep `--features` after wasm-pack's own flags. wasm-pack hands everything from
the first flag it does not recognise to `cargo build`, so a flag or a path
written after `--features wasm` goes to cargo rather than to wasm-pack.

The package is written to `crates/yamc/pkg/`: `yamc.js`, `yamc_bg.wasm` and
their TypeScript declarations. To use it from a page of your own, serve the
directory that holds the page and `pkg/` over HTTP, for instance with
`python -m http.server 8000`, since a JavaScript module does not load from
`file://`.

Core keeps no standalone demo page for this package. The browser page it does
ship is the one `model.to_html(path)` writes (see
[below](#running-a-simulation-in-the-browser-wasmsimulation)), from the
template in `packages/yamc-core/python/yamc/_browser_template/`, with a build
of this crate embedded. In a source build of the Python package,
`scripts/build_wasm_blobs.sh` rebuilds that embedded build (`yamc_sim`) and the
mesh viewer's (`yamt`) into `packages/yamc-core/python/yamc/_wasm/`.

### Running WASM Tests

As core's CI runs them, from `crates/yamc`, with Firefox and geckodriver
installed:

```bash
cd crates/yamc
cargo test -p yamc --features wasm-test
wasm-pack test --headless --firefox --features wasm-test
```

The `wasm-test` feature brings in the WASM-specific tests and their
dependencies. One of them, `tests/wasm_fetch_browser.rs`, fetches nuclear data
from yamc-data.xsplot.com inside the browser, so the run needs network access.


## JavaScript API

The WebAssembly build exposes a JavaScript API that mirrors the Rust API. For example:

```javascript
import init, { Material, Nuclide } from './pkg/yamc.js';

async function main() {
    await init();

    // Create a material
    let material = new Material();
    material.add_element("Li", 1.0);
    material.set_density("g/cm3", 0.534);

    // Calculate the macroscopic total cross section (MT 1) at the
    // default temperature. `macroscopicCrossSection` returns
    // [xs_values, energy_grid]; pass an explicit temperature string
    // (e.g. "300K") as a second argument to override the default.
    let xs_data = material.macroscopicCrossSection(1);
    console.log("Cross section data:", xs_data);
}

main();
```

### Running a simulation in the browser (`WasmSimulation`)

Besides `Material`, the bindings expose a `WasmSimulation` class for
loading a model exported from Python and running the transport loop
in the browser. The exported model JSON (from `model.save(path)`)
carries the geometry, materials, source, and tally definitions; the JS
host supplies the nuclear data files and the run size. (`model.to_html(path)`
bundles that same model into a self-contained, ready-to-run HTML page.)

The exported methods (see the source of truth in
`crates/yamc/src/wasm/simulation_wasm.rs`, and `fetch_wasm.rs` beside it for
`fetchNuclearData`) are:

- `new WasmSimulation()`: construct a simulation; installs a fresh
  in-memory nuclear-data store. Keep a single instance per page.
- `fetchNuclearData(library_url, on_progress)`: fetch every nuclide (and,
  for photon transport, every element) the loaded model needs and the store
  does not hold yet, from the published per-section objects under
  `library_url` (`default_library_url()` names the ENDF/B-VIII.1 set on
  yamc-data.xsplot.com). The library's current release is resolved first
  (`latest.json`, then the manifest it names, checked against it), and every
  section comes from that release. `reactions.arrow` and `energy.arrow` are
  fetched as the byte ranges of the temperature each material is at, or of the
  two that bracket it, so for a material with a temperature set one temperature
  of six comes down the wire; every other section arrives whole and is checked
  against its manifest size and sha256. `on_progress(message, done, total)` is
  optional. Resolves to a JSON summary `{status, library, release,
  manifest_sha256, nuclides, ranged, elements, requests, bytes,
  metadata_bytes}`, where `bytes` counts section data and `metadata_bytes` the
  release pointer and manifest; rejects with the URL that failed, or with the
  loader's own message when a material temperature is outside what the data
  covers.
- `add_file(path, bytes)`: register one virtual nuclear-data file
  (`/<Nuclide>.arrow/<section>`) the host already holds, for instance data
  the exporting Python embedded in the page.
- `file_count()`: number of files currently held in memory.
- `load_model_json(json)`: load a model from its exported JSON,
  replacing any previously loaded model.
- `model_required_nuclides()`: comma-joined list of nuclide names the
  loaded model references.
- `model_missing_nuclides()` / `model_missing_elements()`: the required
  nuclides or photon elements whose section set is not in the store yet,
  comma-joined. Both are empty when transport can run, so a Simulate button
  can be enabled on that.
- `sampleSlice(paramsJson)`: sample a 2D geometry slice; returns a flat
  interleaved `[cell_id, mat_id, ...]` array.
- `boundingBox()`: geometry bounding box as `[cx, cy, cz, wx, wy, wz]`
  (center and widths in cm).
- `geometryJson()`: CsgGeometry JSON for the loaded model, used to build
  the viewer's cell/material legend.
- `plotHtml(paramsJson)`: full interactive geometry-viewer HTML for the
  loaded model.
- `tallyPlotHtml(paramsJson)`: full interactive mesh-tally plot HTML from
  the current tally state.
- `sampleSourcePoints(n, seed)`: sample `n` source positions; returns a
  flat `[x, y, z, ...]` array for the source overlay.
- `simulate_transport(particles, batches, seed)`: run the loaded model
  (overriding its run parameters) and return a JSON result string with a
  per-tally, per-score breakdown of means and standard deviations.
