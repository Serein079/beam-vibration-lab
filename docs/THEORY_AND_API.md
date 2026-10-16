# Theory and implementation

## Model

Two-node cubic Hermite Euler–Bernoulli elements use `[w1, theta1, w2, theta2]`,
where `theta=dw/dx`. The uniform mesh has `n_elements+1` nodes and
`2*(n_elements+1)` global DOFs. SI units are required. The rectangular section
uses `A=b*h` and `I=b*h^3/12`; h is the bending depth.

For element length l:

```text
Ke = EI/l^3 * [[12,6l,-12,6l], [6l,4l²,-6l,2l²],
              [-12,-6l,12,-6l], [6l,2l²,-6l,4l²]]
Me = rhoA*l/420 * [[156,22l,54,-13l], [22l,4l²,13l,-3l²],
                   [54,13l,156,-22l], [-13l,-3l²,-22l,4l²]]
```

`Beam.assemble()` adds these matrices to global indices `[2e,2e+1,2e+2,2e+3]`.
Profiles `EI_profile` [N m²] and `rhoA_profile` [kg/m] override the scalar
properties element by element. Their length, positivity and finiteness are
checked at construction; reconstruct the model when changing its parameters.
Total mass is the sum of `rhoA_e*l`, not simply `rho*A*length` for profiled beams.

## Boundary conditions

| Boundary value | Constrained global DOFs | Physical meaning |
|---|---|---|
| cantilever | 0, 1 | w=theta=0 on left |
| simply-supported | 0, 2n | w=0 on both ends, rotations free |
| clamped-clamped | 0, 1, 2n, 2n+1 | w=theta=0 on both ends |

Homogeneous constraints are eliminated by taking `Kff`, `Mff`, `Ff`.
Nonzero prescribed support displacements are outside this API.
A one-element fixed-fixed model has no free DOFs and is rejected.

## Loads, static displacement and reactions

Uniform q produces consistent element load `[ql/2, ql²/12, ql/2, -ql²/12]`.
`load_vector()` adds a right-end force and dictionaries of node forces/moments.
Only loads placed at mesh nodes are supported; arbitrary physical coordinates
must be represented by a suitable mesh node or independently derived equivalent loads.

`static_result(tip_load=0, distributed_load=0, load=None)` adds scalar loads
and the optional full load vector. It solves `Kff*uf=Ff`. Constrained entries are
zero; support reactions are `(K*u-F)` at constrained DOFs and zero elsewhere.
The result exposes x, displacement, deflection, rotation, reactions and load.
Even reaction entries are forces [N], odd entries are moments [N m].
Applied loads at fixed DOFs appear in reactions.

## Modal analysis

`modal_result(n_modes)` solves `Kff*phi=omega²*Mff*phi` using a symmetric
generalized eigensolver and only requests the lowest needed eigenpairs.
The result's `modes` shape is `(n_modes,n_dof)`, with zero constrained entries.
Normalization is `phi_i.T*M*phi_j=delta_ij`; row signs are chosen so the largest
nodal deflection is positive. Frequencies are omega/(2 pi), in Hz.
Requesting more modes than free DOFs raises an error; nonpositive eigenvalues
are reported rather than silently converted to zero frequency.

Uniform references use `omega=(betaL)²*sqrt(EI/(rhoA*L^4))`:
cantilever roots 1.875104, 4.694091, 7.854757, ...;
simply-supported roots n*pi; fixed-fixed roots 4.730041, 7.853205, 10.995608, ... .
The finite root tables support five modes. Analytical helpers reject profiles,
since uniform formulas would otherwise describe a different model.

## Harmonic response and damping

`harmonic_response(frequencies, load, rayleigh=(alpha_m,alpha_k))` solves:

```text
C = alpha_m*M + alpha_k*K
(K - omega²*M + i*omega*C) U = F
u(t) = Re(U exp(i*omega*t))
```

F may be complex, so multiple force phases can be prescribed. The result has
shape `(n_frequencies,n_dof)`; `abs(U)` is amplitude and `angle(U)` phase.
At zero frequency it reproduces static response. An undamped resonance is
singular or severely ill-conditioned; use positive damping to investigate a peak.
No result clipping or fabricated finite resonance amplitude is applied.

`rayleigh_coefficients([f1,f2],[zeta1,zeta2])` solves two equations:
`zeta_i=alpha_m/(2*omega_i)+alpha_k*omega_i/2`.
Units: alpha_m [1/s], alpha_k [s]. Targets requiring negative coefficients are
rejected because such damping is not nonnegative at all frequencies.

## Transient dynamics

`transient_result()` integrates `M*a+C*v+K*u=F(t)` with constant-step Newmark:

```text
dt=t_final/n_steps
a0=1/(beta*dt²), a1=gamma/(beta*dt)
Keff=K+a0*M+a1*C
```

The effective matrix is factored once. The default beta=1/4, gamma=1/2 is
average acceleration. This API requires gamma>=1/2 and
beta>=(gamma+1/2)²/4, a sufficient unconditional-stability bound for the linear
system. Stability does not guarantee temporal accuracy; refine dt to check
phase and amplitude. Non-default gamma can introduce numerical damping.

Initial displacement comes from released preload or an explicit full vector;
both cannot be supplied together. Initial velocity defaults to zero or may be
provided explicitly. Both vectors must satisfy homogeneous supports.
Initial acceleration is computed from equilibrium at t=0, including `force(0)`.
The force callback must return a finite real global vector of length n_dof.

Results include all displacement/rotation, velocity and acceleration DOFs at
`n_steps+1` times, including the initial state. Energies are
`T=0.5*v.T*M*v`, `V=0.5*u.T*K*u`. Undamped unforced average-acceleration Newmark
conserves T+V up to roundoff. Damped free vibration loses mechanical energy.
Forced vibration can gain energy through external work; decreasing energy is
not a universal validation criterion.

## Compatibility

The original import `import fem_beam as fb` is retained.
`static_deflection()` returns x and nodal w;
`natural_frequencies()` returns Hz and nodal-w mode rows;
`transient_response()` returns time, right-end w and all nodal w.
The right-end value is identically zero for simply-supported or fixed-fixed
beams; use the full result to monitor an interior node.

The dense solver is intended for small meshes. Cubic static tip-load results
can be effectively exact, while very refined eigenproblems can lose accuracy
through conditioning. The linear model omits shear deformation, rotary inertia,
large-deflection effects and 3D stress concentration.
