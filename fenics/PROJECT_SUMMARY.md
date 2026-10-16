# Project Summary

This project models deformation and vibration in beam-like structures with the
finite element method. The simulations use linear elasticity, fixed boundary
constraints, distributed body forces, and surface traction loads.

## Engineering Context

The original work was developed for an aerospace/mechanics project involving
structural vibration and deflection. This cleaned version keeps the core
technical work while removing duplicated assignment files, generated solver
outputs, and temporary build artifacts.

## Methods

- Geometry is represented with FEniCS mesh primitives or CSG geometry.
- Displacement is approximated with first-order vector finite elements.
- Stress is computed from the small-strain linear elastic constitutive model.
- Static problems solve the weak form of equilibrium.
- Dynamic vibration uses a finite-difference time integration scheme initialized
  from a static displacement field.

## Tools

- Python
- FEniCS / DOLFIN
- NumPy
- XDMF/HDF5 output for ParaView visualization

## Interview Talking Points

- How the weak form is assembled from strain energy and external work.
- Why fixed boundary conditions are applied at one end of a cantilever beam.
- How surface traction differs from a point load in finite element modeling.
- How time step size affects stability and resolution in transient vibration.
- Why generated mesh/result files should be excluded from a clean source repo.
