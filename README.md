# Beam Vibration Lab

An educational Python finite element workbench for **Euler–Bernoulli beam bending and vibration**.
Analyze supports, static loads, natural frequencies, forced frequency response,
and transient energy using NumPy and SciPy. The original cantilever API remains available.

[中文说明](README.zh-CN.md) · [Theory and API](docs/THEORY_AND_API.md) · [Validation](docs/VALIDATION.md) · [Publishing](docs/PUBLISHING.zh-CN.md)

**Source:** an extension of [yanyana117/finite-element-vibration-simulation](https://github.com/yanyana117/finite-element-vibration-simulation),
baseline `8eb430f`. The upstream license is unspecified; read [NOTICE.md](NOTICE.md) before public redistribution.

| Support comparison | Frequency response | Mechanical energy |
|---|---|---|
| ![supports](results/support_comparison.png) | ![response](results/harmonic_response.png) | ![energy](results/energy_check.png) |

## What this edition adds

- Cantilever, simply-supported and clamped-clamped supports.
- Nodal forces and moments; full static deflections, rotations and support reactions.
- Mass-normalized full modal vectors and analytical references for three supports.
- Complex harmonic response with Rayleigh damping; two-frequency damping calibration.
- Newmark displacement, velocity and acceleration histories; kinetic and strain energies.
- JSON case files and `beam-fem` CLI exporting CSV, NPZ and a replayable input file.
- Installable package, numerical regression tests and Windows/Linux CI configuration.

The mesh is uniform. Elementwise EI and rhoA profiles model stepped sections.
The 3D legacy FEniCS scripts are optional and have **not** been validated in the current environment.

## Install and run

Requires Python 3.11 or newer. From the repository directory:

```bash
python -m venv .venv
# Windows PowerShell: .\.venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
beam-fem --config configs/cantilever_static.json --output runs/static
python examples/extended_studies.py
python scripts/run_all.py --output runs/latest
```

Core-only installation: `python -m pip install .`. Figures require `.[plot]`.
The full workflow runs all tests, seven example scripts and four CLI cases,
and writes logs, JUnit XML, a manifest and Markdown report. A failure produces a nonzero exit code.

## Python API

```python
from fem_beam import steel_rectangular, rayleigh_coefficients

beam = steel_rectangular(n_elements=30, boundary="simply-supported")
load = beam.load_vector(nodal_forces={15: -1.0})
static = beam.static_result(load=load)
print(static.deflection[15], static.reactions)

modal = beam.modal_result(4)
damping = rayleigh_coefficients(modal.frequencies[:2], [0.02, 0.02])
response = beam.harmonic_response([0, 20, 30], load, rayleigh=damping)
transient = beam.transient_result(0.5, 2000,
    initial_displacement=static.displacement, rayleigh=damping)
print(transient.total_energy[-1] / transient.total_energy[0])
```

All inputs and arrays use SI units. DOFs alternate `[w0, theta0, w1, theta1, ...]`.
Positive rotation is `dw/dx`. The complex harmonic convention is `Re(U exp(i omega t))`.
CLI transient cases release the specified static preload at `t=0`; the Python API also accepts `force(t)`.

## Repository map

```text
src/          solver and CLI
configs/      four replayable analysis cases
examples/     six upstream studies and one extended study
tests/        analytical checks, dynamic invariants and CLI checks
scripts/      one-command validation workflow
results/      reproducible example figures and extension data
docs/         theory, verification and publishing instructions
fenics/       optional legacy 3D studies (separate environment)
.github/      CI and issue / pull request templates
```

## Scope

Small-displacement linear elasticity, slender beams, no shear deformation,
rotary inertia, plasticity or contact. Dense matrices target small educational
models rather than large production meshes. Very fine meshes can amplify roundoff;
frequency convergence is not guaranteed to remain monotonic. The original material
study includes large deflections and is a linear-theory demonstration, not a nonlinear design prediction.

The checked-in validation record describes the actual local Windows/Python run.
CI is configured for Python 3.11–3.13 on Windows and Linux; remote CI results require
pushing to GitHub. See [CONTRIBUTING.md](CONTRIBUTING.md) for adding new models.
