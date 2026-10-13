# Installation

YAMC can be installed and used on Linux, Windows and MacOS.

## Standard Installation (Recommended)

YAMC is a Rust back end package with a Python API. Installation requires no additional Python package dependencies.
YAMC is available on [PyPI](https://pypi.org/project/yamc).

### Install from PyPI

To install, make sure you have [pip installed](https://pip.pypa.io/en/stable/installation/), then run:

```bash
pip install yamc
```

This installs the core YAMC transport solver with multi-threading support (via Rayon). MPI is not required.



### Optional extras

YAMC provides optional extras for additional functionality:

| Extra | Install command | What it adds |
|-------|----------------|--------------|
| cad | `pip install yamc[cad]` | CAD meshing pipeline, surface and tetrahedral (requires cadquery) |
| viz | `pip install yamc[viz]` | VTKHDF visualization export (h5py, numpy) |
| pandas | `pip install yamc[pandas]` | pandas DataFrame export of tally results |
| polars | `pip install yamc[polars]` | polars DataFrame export of tally results |
| xarray | `pip install yamc[xarray]` | xarray (labelled N-D array) export of tally results |
| all | `pip install yamc[all]` | cad and viz (does not include pandas, polars or xarray) |
| dev | `pip install yamc[dev]` | Development dependencies (pytest, maturin, cadquery) |

To install everything, including the CAD-to-mesh-to-transport pipeline:

```bash
pip install yamc[all]
```

### Virtual Environment (Recommended)

Further instructions are in the [Python docs](https://packaging.python.org/en/latest/guides/installing-using-pip-and-virtual-environments/).

```bash
python3 -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
pip install yamc[all]
```

### Install the development version

To get the latest unreleased changes, install directly from the GitHub repository. This builds from source, so it needs a Rust toolchain, which can be installed with [rustup](https://www.rust-lang.org/tools/install):

```bash
# Install Rust (skip if you already have it)
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

pip install git+https://github.com/fusion-neutronics/yamc.git
```

Extras work too, e.g. `pip install "yamc[all] @ git+https://github.com/fusion-neutronics/yamc.git"`.


## MPI Installation (Advanced)

YAMC can be built with MPI support for distributed-memory parallelism across multiple nodes, for example on large clusters.

MPI is optional, and the standard installation above is enough for most users. Install MPI only if you need:
- Multi-node parallelism (HPC clusters)
- Scaling beyond a single machine
- Hybrid MPI + threading parallelism

### Prerequisites

First, install an MPI implementation on your system:

#### Ubuntu/Debian

```bash
sudo apt-get update && sudo apt-get install -y \
    libopenmpi-dev openmpi-bin \
    build-essential clang libclang-dev
```

#### Fedora/RHEL/CentOS

```bash
sudo dnf install openmpi-devel
# or
sudo yum install openmpi-devel

# Load MPI module (if using environment modules)
module load mpi/openmpi-x86_64
```

#### macOS (Homebrew)

```bash
brew install open-mpi
```

#### Windows (Microsoft MPI)

1. Download and install [Microsoft MPI](https://learn.microsoft.com/en-us/message-passing-interface/microsoft-mpi)
2. Install both `msmpisetup.exe` (runtime) and `msmpisdk.msi` (SDK)

#### Verify MPI Installation

```bash
# Check MPI compiler is available
mpicc --version

# Test MPI works
mpirun -np 2 hostname
```

### Install YAMC with MPI

Building with MPI requires the Rust toolchain. To install Rust:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
source $HOME/.cargo/env
```

Then clone and build from source with the `mpi` feature:

```bash
git clone https://github.com/fusion-neutronics/yamc.git
cd yamc
pip install maturin
maturin develop --release --features mpi
```

### Verify MPI Installation

```bash
# Check MPI is available in YAMC
python -c "import yamc; print(f'MPI enabled: {yamc.parallel.has_mpi()}')"
# Should print: MPI enabled: True

# Test with multiple ranks
mpirun -np 2 python -c "import yamc; print(f'Rank {yamc.parallel.mpi_rank()} of {yamc.parallel.mpi_size()}')"
# Should print:
# Rank 0 of 2
# Rank 1 of 2
```

## Docker Installation

### Standard Docker Image

```dockerfile
FROM python:3.11-slim

RUN pip install yamc[all]

WORKDIR /workspace
CMD ["python"]
```

### Docker Image with MPI

```dockerfile
FROM ubuntu:22.04

# Install system dependencies
RUN apt-get update && apt-get install -y \
    python3 python3-pip \
    libopenmpi-dev openmpi-bin \
    build-essential clang libclang-dev curl git \
    && rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Install YAMC with MPI from source
RUN git clone https://github.com/fusion-neutronics/yamc.git /tmp/yamc \
    && pip3 install maturin \
    && cd /tmp/yamc \
    && maturin develop --release --features mpi \
    && rm -rf /tmp/yamc

WORKDIR /workspace
CMD ["python3"]
```

Build and run:

```bash
docker build -t yamc-mpi -f Dockerfile .
docker run --rm yamc-mpi mpirun -np 4 python3 -c "import yamc; print(yamc.parallel.mpi_rank())"
```

## HPC/Cluster Installation

For HPC environments with environment modules:

```bash
# Load required modules
module load python/3.11
module load openmpi/4.1.4
module load rust/1.75.0  # or install Rust locally

# Create virtual environment
python3 -m venv ~/.venvs/yamc-mpi
source ~/.venvs/yamc-mpi/bin/activate

# Clone and build with MPI
git clone https://github.com/fusion-neutronics/yamc.git
cd yamc
pip install maturin
maturin develop --release --features mpi

# Test installation
mpirun -np 4 python -c "import yamc; print(f'Rank {yamc.parallel.mpi_rank()}')"
```

See [simulating](simulate_transport.md#running-the-simulation) for MPI usage, hybrid MPI+threading, and Slurm examples.
