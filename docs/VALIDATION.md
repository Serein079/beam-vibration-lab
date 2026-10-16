# Validation method

The tests verify physical quantities and invariants, rather than comparing the
implementation to a duplicate of its own algorithm. Actual run evidence is in
[LOCAL_VALIDATION.md](LOCAL_VALIDATION.md), generated after the workflow completes.

| Check | Independent reference |
|---|---|
| Cantilever end force | PL³/(3EI) |
| Cantilever uniform load | qL⁴/(8EI) |
| Cantilever end moment | w=ML²/(2EI), theta=ML/EI |
| Simply-supported uniform load | w(mid)=5qL⁴/(384EI), reactions=-qL/2 |
| Simply-supported midpoint force | w(mid)=PL³/(48EI) |
| Fixed-fixed uniform load | w(mid)=qL⁴/(384EI), end moments ±qL²/12 |
| Three support types, first four frequencies | Euler–Bernoulli characteristic roots |
| Modal vectors | eigen-equilibrium residual; Phi*M*Phi.T=identity |
| Harmonic zero frequency | static solution |
| Complex harmonic response | global dynamic-equilibrium residual |
| Modal harmonic force at resonance | single-oscillator response with known damping |
| Newmark free, undamped | conserved mechanical energy |
| Newmark free, damped | nonincreasing mechanical energy |
| Initial modal velocity | sinusoidal single-oscillator solution |
| Constant force and matching static initial condition | unchanged equilibrium |
| Rayleigh fit | independently recover the two target damping ratios |
| CLI | exports, input replay, exit status and typo rejection |

Invalid geometry, profile lengths, material values, node indices, mode counts,
force-vector sizes and time-integration parameters are checked explicitly.
Nine upstream tests remain and are supplemented by extension and CLI tests.

Run `python -m pytest -q`, or run `python scripts/run_all.py --output runs/latest`
for tests, seven example scripts and four JSON analyses. Logs, test XML and
exit codes are preserved. `results/extended_summary.json` records support
comparisons and energy quantities for a small-deflection -1 N preload.

The historical upstream material study uses -50 N on a 0.5 m beam, leading to
deflection/length ratios above 0.1. Those outputs illustrate linear parameter
scaling and cannot establish real nonlinear structural behavior. The new
demonstrations use smaller loads.

The original cut-out example has four weakened regions (12 of 60 elements):
local EI is multiplied by 0.4 and rhoA by 0.6. The whole-beam mass ratio is
0.92. It represents an equivalent section profile, not explicit holes.

Frequency mesh convergence should be inspected over its resolved range.
Roundoff can reverse error reduction at very fine discretization; this edition
does not claim universal monotonic fourth-order convergence.

**Legacy FEniCS:** Python syntax is checked, but no 3D mesh or solution is
generated in the local validation run. Fixes to the box geometry, body-force
units, traction direction and transient timestamps remain runtime-unverified.
The implicit legacy time scheme is distinct from the 1D Newmark implementation.
