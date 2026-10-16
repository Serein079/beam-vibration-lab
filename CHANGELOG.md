# Changelog

## 0.2.0 — local extension edition, 2026-10-09

- Generalize homogeneous supports to cantilever, simply-supported and fixed-fixed.
- Add construction/input validation and physically meaningful error messages.
- Add node forces/moments, full static rotations and support reactions, total mass.
- Add full mass-normalized modes and selected eigenpair computation.
- Add uniform analytical modal references for simply-supported and fixed-fixed beams.
- Add complex harmonic response and two-frequency Rayleigh damping calibration.
- Add initial-condition controls, full dynamic DOFs and kinetic/strain energy histories.
- Retain the original three tuple-returning APIs and all six original studies.
- Add four JSON cases, CLI exports/replay, packaging, workflow and Windows/Linux CI.
- Add analytical, modal, harmonic, transient and CLI validation tests.
- Correct 5x10 mm captions to the modeled 50x10 mm section.
- Correct the weakened-section description to four regions, 60% local EI reduction,
  40% local rhoA reduction and 8% whole-beam mass reduction.
- Legacy FEniCS: change the mislabeled side-void section to a true closed box;
  use rho*g body force; define top traction downward in z; write t=0 initial
  state; initialize the prior state using release acceleration. Runtime unverified.
- Add English/Chinese usage, theory, scope, validation, contribution and publishing docs.
- Record upstream provenance without inventing an upstream license.

## Baseline

Upstream commit `8eb430f666c048604d37296c82a0364c26fe4117`:
cantilever statics, modes and Newmark response, material/profile studies,
nine validation tests and three legacy FEniCS examples.
