# Contributing

Install development dependencies with `python -m pip install -e '.[dev]'`.
Run `python scripts/run_all.py --output runs/development` before submitting changes.
Build with `python -m build --sdist --wheel`. Keep generated caches and virtual environments out of Git.

Each new physical feature should include an independent analytical reference,
an equilibrium/conservation invariant, or a trusted external benchmark.
Document units, sign conventions, boundary conditions and applicable assumptions.
Do not lower tolerances simply to hide a failure; explain conditioning or discretization limits.
Preserve original APIs where possible and update the changelog for deliberate changes.

Submit focused pull requests with the physical problem, resulting behavior and
actual validation evidence. Legacy FEniCS changes require tests in the separate
backend before they can be described as runtime-validated.
Read NOTICE.md: the upstream license has not been specified, and this project
does not resolve redistribution rights by assigning a blanket license.
