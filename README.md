# yamc

**Y**et **A**nother **M**onte **C**arlo: an open-source Monte Carlo particle
transport code for neutron and photon simulations, with a focus on fusion
neutronics workflows.

[Documentation](https://fusion-neutronics.github.io/yamc/) |
[PyPI](https://pypi.org/project/yamc/)

```bash
pip install yamc
```

```python
import yamc

yamc.cross_section_data = "endf-b8.1"

steel = yamc.Material({"Fe56": 1.0}, density=7.87)
sphere = yamc.Sphere(r=10.0)
cell = yamc.Cell(region=-sphere, fill=steel)
geometry = yamc.Geometry([cell])

source = yamc.NeutronSource(energy=yamc.Discrete([14.06e6], [1.0]))
tally = yamc.Tally(scores=["flux"])

model = yamc.Model(geometry=geometry, source=source, tallies=[tally])
results = model.simulate_transport(total_particles=10_000)
print(results[tally].mean)
```

## What it does

Continuous-energy neutron and photon transport over CSG or unstructured mesh
geometry, on CPU or GPU, with tallies, variance reduction, transport-coupled
transmutation and shutdown dose rates. The solver is Rust; the wheel needs no
compiler and no external toolchain.

## What this repository contains

Documentation, packaging and release automation. The code lives in the
`yamc-core` wheel this package depends on, which provides the `yamc` module
itself, its type stubs and its docstrings. `pip install yamc` pulls it in.

That arrangement is why there is no Python here: one source of truth for the
API, and no re-export layer that can drift from it.

## Transmutation without transport

If you already have a spectrum and want none of the transport stack, install
[yani](https://github.com/fusion-neutronics/yani) instead: same materials, same
schedules, same nuclear data, no geometry.

The two are **alternatives, not companions**. Installing both puts two extension
modules in one process, so `yamc.Material` and `yani.Material` are distinct
types and the nuclear-data configuration exists twice.

## License

MIT.
