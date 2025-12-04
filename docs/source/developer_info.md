# Developer Info

This guide shows how to clone, build, test, and generate documentation for the project. Instructions are written for Linux but equivilent commands should be possible in Windows and Mac OS.

## 1. Clone the Repository
```bash
git clone https://github.com/fusion-neutronics/materials_for_mc.git
cd materials_for_mc
```

## 2. Build & Test the Rust Library
Install Rust Linux and Mac: https://www.rust-lang.org/tools/install
```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
```
Windows https://forge.rust-lang.org/infra/other-installation-methods.html

Install the nessecary build tools
```
sudo apt update && sudo apt install build-essential libssl-dev pkg-config
```

Build Rust Package:
```bash
cargo build
```
Run Rust tests:
```bash
cargo test
```

## 3. Build the Python Extension (maturin)

(Optional) Create a Python virtual environment:
```bash
sudo apt install python3-venv
python3 -m venv .venv
source .venv/bin/activate  # Windows: .venv\Scripts\activate
```

Install maturin (once):
```bash
pip install maturin
```

Build & develop-install the Python module in debug mode:
```bash
maturin develop --features pyo3
```
Build & develop-install the Python module in release mode for a significant speed improvement
```bash
maturin develop --release --features pyo3
```

(Optional) Build a wheel:
```bash
maturin build --features pyo3
```

## 4. Run Python Tests
Install test dependencies:
```bash
pip install pytest
```
Run tests:
```bash
pytest
```

## 5. Build Rust API Docs
```bash
cargo doc --no-deps
# Open target/doc/materials_for_mc/index.html in a browser
```

## 6. Build Python (Sphinx) Docs
Install doc requirements:
```bash
pip install -r docs/requirements.txt
```
Ensure extension is built (maturin develop) then build docs:
```bash
sphinx-build -b html docs/source docs/build/html
```
Open docs:
```bash
python -m webbrowser docs/build/html/index.html
```


### Profile the code

High level profile
```
# Build with debug symbols in release mode
cargo build --release --example tbr
RUSTFLAGS="-C force-frame-pointers=yes" cargo build --release --example tbr

# Profile with full call graph
sudo perf record -g --call-graph dwarf ./target/release/examples/tbr

# View interactive report
sudo perf report --no-children

# Or generate text report with full call stacks
sudo perf report --no-children --stdio > perf_report.txt
```

Detailed profile
```
RUSTFLAGS="-C force-frame-pointers=yes" cargo build --release --example tbr
sudo perf record -g --call-graph dwarf ./target/release/examples/tbr
sudo perf report --no-children --stdio | head -200
```