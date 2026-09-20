# WebAssembly

The yamc library provides WebAssembly (WASM) bindings, allowing the library to run in web browsers and other JavaScript environments. This enables client-side nuclear data processing and visualization without requiring server-side computation.

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

### Prerequisites

First, install the required tools:

```bash
curl https://rustwasm.github.io/wasm-pack/installer/init.sh -sSf | sh
```

### Building the Package

To quickly build the WebAssembly package:

```bash
wasm-pack build --target web --features wasm crates/yamc
```

To build the optimized WebAssembly package:
```bash
wasm-pack build --target web --features wasm --release crates/yamc
```

This creates a `pkg/` directory containing the generated WebAssembly files and JavaScript bindings.

Create an index.html that makes use of the WASM file.

Copy the generated package to the examples directory and serve locally:

```bash
cp -r pkg examples/wasm/
# Serve the demo pages locally
python -m http.server 8000
# Open the demo pages in your browser
firefox http://localhost:8000/examples/wasm/index.html
```

### Running WASM Tests

To test the WebAssembly bindings:

```bash
wasm-pack test --headless --firefox crates/yamc
```

### Running WASM-Specific Tests

For tests that require WASM-specific features, enable the optional `wasm-test` feature:

```bash
cargo test -p yamc --features wasm-test
wasm-pack test --headless --firefox --features wasm-test crates/yamc
```

This ensures that WASM-specific tests and dependencies are only included when needed.


## JavaScript API

The WebAssembly build exposes a JavaScript API that mirrors the Rust API. Example usage:

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

Beyond `Material`, the bindings expose a `WasmSimulation` class for
loading a model exported from Python and running the transport loop
in the browser. The exported model JSON (from `model.save(path)`)
carries the geometry, materials, source, and tally definitions; the JS
host supplies the nuclear data files and the run size. (`model.to_html(path)`
bundles that same model into a self-contained, ready-to-run HTML page.)

The exported methods (see the source of truth in
`crates/yamc/src/wasm/simulation_wasm.rs`) are:

- `new WasmSimulation()` -- construct a simulation; installs a fresh
  in-memory nuclear-data store. Keep a single instance per page.
- `fetchNuclearData(library_url, on_progress)` -- fetch every nuclide (and,
  for photon transport, every element) the loaded model needs and the store
  does not hold yet, from the published per-section objects under
  `library_url` (`default_library_url()` names the ENDF/B-VIII.1 set on
  yamc-data.xsplot.com). `reactions.arrow` and `energy.arrow` are fetched as
  the byte ranges of the temperature each material is at, or of the two that
  bracket it, so one temperature of six comes down the wire; every other
  section arrives whole. `on_progress(message, done, total)` is optional.
  Resolves to a JSON summary `{status, nuclides, ranged, elements, requests,
  bytes}`; rejects with the URL that failed, or with the loader's own message
  when a material temperature is outside what the data covers.
- `add_file(path, bytes)` -- register one virtual nuclear-data file
  (`/<Nuclide>.arrow/<section>`) the host already holds, for instance data
  the exporting Python embedded in the page.
- `file_count()` -- number of files currently held in memory.
- `load_model_json(json)` -- load a model from its exported JSON,
  replacing any previously loaded model.
- `model_required_nuclides()` -- comma-joined list of nuclide names the
  loaded model references.
- `model_missing_nuclides()` / `model_missing_elements()` -- the required
  nuclides / photon elements whose section set is not in the store yet,
  comma-joined. Empty when transport can run; what a Simulate button gates on.
- `sampleSlice(paramsJson)` -- sample a 2D geometry slice; returns a flat
  interleaved `[cell_id, mat_id, ...]` array.
- `boundingBox()` -- geometry bounding box as `[cx, cy, cz, wx, wy, wz]`
  (center + widths in cm).
- `geometryJson()` -- CsgGeometry JSON for the loaded model, used to build
  the viewer's cell/material legend.
- `plotHtml(paramsJson)` -- full interactive geometry-viewer HTML for the
  loaded model.
- `tallyPlotHtml(paramsJson)` -- full interactive mesh-tally plot HTML from
  the current tally state.
- `sampleSourcePoints(n, seed)` -- sample `n` source positions; returns a
  flat `[x, y, z, ...]` array for the source overlay.
- `simulate_transport(particles, batches, seed)` -- run the loaded model
  (overriding its run parameters) and return a JSON result string with a
  per-tally, per-score breakdown of means and standard deviations.
