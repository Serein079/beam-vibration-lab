"""
fem_beam
========

A finite element solver for Euler-Bernoulli beam bending and vibration,
implemented in NumPy/SciPy.

The beam is discretised with 2-node cubic (Hermite) elements. Each node carries
two degrees of freedom: transverse deflection ``w`` and rotation ``theta``.
The element stiffness and consistent-mass matrices are the standard closed-form
4x4 matrices derived from the cubic shape functions:

    Ke = (EI / L^3) *  [[ 12,   6L,  -12,   6L ],
                        [ 6L, 4L^2,  -6L, 2L^2 ],
                        [-12,  -6L,   12,  -6L ],
                        [ 6L, 2L^2,  -6L, 4L^2 ]]

    Me = (rhoA L / 420) * [[ 156,   22L,    54,  -13L ],
                           [ 22L,  4L^2,   13L, -3L^2 ],
                           [  54,   13L,   156,  -22L ],
                           [-13L, -3L^2,  -22L,  4L^2 ]]

Supported analyses:
  * static bending          K u = f
  * modal / free vibration  K phi = omega^2 M phi
  * transient dynamics      M a + C v + K u = f(t)   (Newmark-beta)

The section can be uniform or vary element-by-element (via ``EI_profile`` /
``rhoA_profile``), which captures stepped or cut-out beams. Validated against
analytical solutions in tests/.
"""
from dataclasses import dataclass
from typing import Optional, Callable

import numpy as np
from scipy.linalg import eigh

__version__ = "0.2.0"


@dataclass
class StaticResult:
    """Full displacement/rotation vector and constrained-DOF reactions (SI)."""
    x: np.ndarray
    displacement: np.ndarray
    reactions: np.ndarray
    load: np.ndarray

    @property
    def deflection(self):
        return self.displacement[0::2]

    @property
    def rotation(self):
        return self.displacement[1::2]


@dataclass
class ModalResult:
    """Frequencies [Hz] and full, mass-normalized modes (one mode per row)."""
    frequencies: np.ndarray
    modes: np.ndarray


@dataclass
class TransientResult:
    """Full DOF histories, mechanical energies [J], and time [s]."""
    time: np.ndarray
    displacement: np.ndarray
    velocity: np.ndarray
    acceleration: np.ndarray
    kinetic_energy: np.ndarray
    strain_energy: np.ndarray

    @property
    def deflection(self):
        return self.displacement[:, 0::2]

    @property
    def total_energy(self):
        return self.kinetic_energy + self.strain_energy


def _positive(value, name):
    if not np.isscalar(value) or not np.isfinite(value) or value <= 0:
        raise ValueError(f"{name} must be finite and positive")


def _count(value, name):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)) or value < 1:
        raise ValueError(f"{name} must be a positive integer")

# Dimensionless eigenvalues (beta_n * L) for a clamped-free Euler-Bernoulli beam
CANTILEVER_BETA_L = np.array([1.8751040687, 4.6940911330, 7.8547574382,
                              10.9955407349, 14.1371683910])


@dataclass
class Beam:
    """An Euler-Bernoulli beam with uniform mesh and selectable supports.

    Parameters
    ----------
    length : float        beam length L [m]
    E : float             Young's modulus [Pa]
    I : float             second moment of area [m^4]
    rho : float           density [kg/m^3]
    A : float             cross-sectional area [m^2]
    n_elements : int      number of finite elements
    EI_profile : array    optional per-element bending stiffness E*I (len n_elements)
    rhoA_profile : array  optional per-element mass-per-length rho*A (len n_elements)
    boundary : str       cantilever, simply-supported, or clamped-clamped
    """
    length: float
    E: float
    I: float
    rho: float
    A: float
    n_elements: int = 20
    EI_profile: Optional[np.ndarray] = None
    rhoA_profile: Optional[np.ndarray] = None
    boundary: str = "cantilever"

    def __post_init__(self):
        for name in ("length", "E", "I", "rho", "A"):
            _positive(getattr(self, name), name)
        _count(self.n_elements, "n_elements")
        if self.boundary not in ("cantilever", "simply-supported", "clamped-clamped"):
            raise ValueError("boundary must be cantilever, simply-supported or clamped-clamped")
        if self.boundary == "clamped-clamped" and self.n_elements < 2:
            raise ValueError("clamped-clamped requires at least two elements")
        for name in ("EI_profile", "rhoA_profile"):
            profile = getattr(self, name)
            if profile is not None:
                profile = np.asarray(profile, dtype=float)
                if profile.shape != (self.n_elements,) or not np.all(np.isfinite(profile)) or np.any(profile <= 0):
                    raise ValueError(f"{name} must contain n_elements finite positive values")
                setattr(self, name, profile.copy())

    @property
    def constrained_dofs(self):
        if self.boundary == "cantilever":
            return np.array([0, 1])
        if self.boundary == "simply-supported":
            return np.array([0, self.n_dof - 2])
        return np.array([0, 1, self.n_dof - 2, self.n_dof - 1])

    @property
    def total_mass(self):
        return sum(self._rhoA_e(e) * self.le for e in range(self.n_elements))

    def _vector(self, value, name, complex_allowed=False):
        if not complex_allowed and np.iscomplexobj(value):
            raise ValueError(f"{name} must be real-valued")
        vector = np.asarray(value, dtype=complex if complex_allowed else float)
        if vector.shape != (self.n_dof,) or not np.all(np.isfinite(vector)):
            raise ValueError(f"{name} must be a finite vector of length {self.n_dof}")
        return vector

    # --- geometry helpers ---
    @property
    def n_nodes(self):
        return self.n_elements + 1

    @property
    def n_dof(self):
        return 2 * self.n_nodes

    @property
    def le(self):
        return self.length / self.n_elements

    @property
    def node_x(self):
        return np.linspace(0.0, self.length, self.n_nodes)

    def _EI_e(self, e):
        if self.EI_profile is not None:
            return float(self.EI_profile[e])
        return self.E * self.I

    def _rhoA_e(self, e):
        if self.rhoA_profile is not None:
            return float(self.rhoA_profile[e])
        return self.rho * self.A

    # --- element matrices ---
    def _element_stiffness(self, EI):
        L = self.le
        return (EI / L ** 3) * np.array([
            [12,    6 * L,   -12,    6 * L],
            [6 * L, 4 * L * L, -6 * L, 2 * L * L],
            [-12,  -6 * L,    12,   -6 * L],
            [6 * L, 2 * L * L, -6 * L, 4 * L * L],
        ])

    def _element_mass(self, rhoA):
        L = self.le
        return (rhoA * L / 420.0) * np.array([
            [156,    22 * L,    54,   -13 * L],
            [22 * L, 4 * L * L, 13 * L, -3 * L * L],
            [54,     13 * L,   156,   -22 * L],
            [-13 * L, -3 * L * L, -22 * L, 4 * L * L],
        ])

    # --- assembly ---
    def assemble(self):
        """Assemble the global stiffness K and consistent mass M matrices."""
        n = self.n_dof
        K = np.zeros((n, n))
        M = np.zeros((n, n))
        for e in range(self.n_elements):
            dofs = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
            K[np.ix_(dofs, dofs)] += self._element_stiffness(self._EI_e(e))
            M[np.ix_(dofs, dofs)] += self._element_mass(self._rhoA_e(e))
        return K, M

    def _free_dofs(self):
        return np.setdiff1d(np.arange(self.n_dof), self.constrained_dofs)

    def load_vector(self, tip_load=0.0, distributed_load=0.0,
                    nodal_forces=None, nodal_moments=None):
        """Equivalent nodal loads. Forces [N], moments [N m], uniform q [N/m].

        ``nodal_forces`` and ``nodal_moments`` map zero-based node indices to
        loads. Positive displacement is positive force; theta = dw/dx.
        Loading a constrained DOF changes its reaction, not its displacement.
        """
        if not np.all(np.isfinite([tip_load, distributed_load])):
            raise ValueError("loads must be finite")
        f = self._load_vector(tip_load, distributed_load)
        for mapping, offset in ((nodal_forces, 0), (nodal_moments, 1)):
            for node, value in (mapping or {}).items():
                if isinstance(node, bool) or not isinstance(node, (int, np.integer)) or not 0 <= node < self.n_nodes:
                    raise ValueError("load node must be an integer in [0, n_nodes)")
                if not np.isfinite(value):
                    raise ValueError("nodal load must be finite")
                f[2 * node + offset] += value
        return f

    def static_result(self, tip_load=0.0, distributed_load=0.0, load=None):
        """Solve K u = f; optional full global ``load`` is added to scalar loads."""
        K, _ = self.assemble()
        f = self.load_vector(tip_load, distributed_load)
        if load is not None:
            f += self._vector(load, "load")
        free = self._free_dofs()
        u = np.zeros(self.n_dof)
        u[free] = np.linalg.solve(K[np.ix_(free, free)], f[free])
        reactions = np.zeros(self.n_dof)
        reactions[self.constrained_dofs] = (K @ u - f)[self.constrained_dofs]
        return StaticResult(self.node_x, u, reactions, f)

    def _load_vector(self, tip_load=0.0, distributed_load=0.0):
        f = np.zeros(self.n_dof)
        if distributed_load:
            L = self.le
            q = distributed_load
            fe = np.array([q * L / 2, q * L * L / 12, q * L / 2, -q * L * L / 12])
            for e in range(self.n_elements):
                dofs = [2 * e, 2 * e + 1, 2 * e + 2, 2 * e + 3]
                f[dofs] += fe
        if tip_load:
            f[2 * (self.n_nodes - 1)] += tip_load
        return f

    # --- static analysis ---
    def _static_full(self, tip_load=0.0, distributed_load=0.0):
        return self.static_result(tip_load, distributed_load).displacement

    def static_deflection(self, tip_load=0.0, distributed_load=0.0):
        """Solve static bending; return node x-coordinates and deflection w."""
        u = self._static_full(tip_load, distributed_load)
        return self.node_x, u[0::2]

    # --- modal analysis ---
    def modal_result(self, n_modes=4):
        """Lowest modes with full DOFs and M-normalization; errors on excess modes."""
        _count(n_modes, "n_modes")
        K, M = self.assemble()
        free = self._free_dofs()
        if n_modes > len(free):
            raise ValueError(f"n_modes must not exceed {len(free)} free DOFs")
        eigvals, eigvecs = eigh(K[np.ix_(free, free)], M[np.ix_(free, free)],
                                subset_by_index=(0, n_modes - 1))
        if np.any(eigvals <= 0):
            raise ValueError("nonpositive eigenvalue: check model conditioning")
        freqs = np.sqrt(eigvals) / (2 * np.pi)
        modes = []
        for k in range(n_modes):
            full = np.zeros(self.n_dof)
            full[free] = eigvecs[:, k]
            shape = full[0::2]
            if shape[np.argmax(np.abs(shape))] < 0:
                full = -full
            modes.append(full)
        return ModalResult(freqs, np.array(modes))

    def natural_frequencies(self, n_modes=4):
        """Compatibility API: frequencies [Hz], nodal deflections of mode shapes."""
        result = self.modal_result(n_modes)
        return result.frequencies, result.modes[:, 0::2]

    def harmonic_response(self, frequencies, load, rayleigh=(0.0, 0.0)):
        """Complex displacement phasors for f(t)=Re(F exp(i omega t)).

        Returns (n_frequencies, n_dof). Rayleigh damping regularizes resonance;
        undamped response exactly at an eigenfrequency is singular.
        """
        frequencies = np.asarray(frequencies, dtype=float)
        if frequencies.ndim != 1 or len(frequencies) == 0 or not np.all(np.isfinite(frequencies)) or np.any(frequencies < 0):
            raise ValueError("frequencies must be a nonempty vector of finite nonnegative Hz")
        am, ak = self._damping(rayleigh)
        f = self._vector(load, "load", complex_allowed=True)
        K, M = self.assemble()
        free = self._free_dofs()
        Kr, Mr = K[np.ix_(free, free)], M[np.ix_(free, free)]
        Cr = am * Mr + ak * Kr
        response = np.zeros((len(frequencies), self.n_dof), dtype=complex)
        for i, frequency in enumerate(frequencies):
            omega = 2 * np.pi * frequency
            response[i, free] = np.linalg.solve(Kr - omega**2 * Mr + 1j * omega * Cr, f[free])
        return response

    @staticmethod
    def _damping(rayleigh):
        coefficients = np.asarray(rayleigh, dtype=float)
        if coefficients.shape != (2,) or not np.all(np.isfinite(coefficients)) or np.any(coefficients < 0):
            raise ValueError("rayleigh must contain two finite nonnegative coefficients")
        return coefficients

    # --- transient dynamics (Newmark-beta) ---
    def transient_result(self, t_final, n_steps, preload_tip=0.0,
                           preload_distributed=0.0, force=None,
                           rayleigh=(0.0, 0.0), beta=0.25, gamma=0.5,
                           initial_displacement=None, initial_velocity=None):
        """Integrate  M a + C v + K u = f(t)  with the Newmark-beta method.

        Parameters
        ----------
        t_final : float            total simulated time [s]
        n_steps : int              number of time steps
        preload_tip,
        preload_distributed        static load whose deflection is used as the
                                   initial condition u0 (then released unless a
                                   ``force`` keeps acting). Models a beam pulled
                                   to a static shape and let go (ring-down).
        force : callable(t)->vec   optional time-varying global force vector
                                   (length n_dof). Default: free vibration.
        rayleigh : (alpha_m, alpha_k)  Rayleigh damping  C = am*M + ak*K
        beta, gamma : float        Newmark parameters (default = average
                                   acceleration: unconditionally stable).

        Returns
        -------
        TransientResult
            Full DOF histories (n_steps+1, n_dof), including rotations, and
            kinetic/strain energies. Explicit initial conditions must satisfy
            supports. ``initial_displacement`` and preload are alternatives.
        """
        _positive(t_final, "t_final")
        _count(n_steps, "n_steps")
        _positive(beta, "beta")
        _positive(gamma, "gamma")
        if gamma < 0.5 or beta < 0.25 * (gamma + 0.5)**2:
            raise ValueError("require gamma >= 0.5 and beta >= (gamma + 0.5)^2 / 4 for unconditional stability")
        K, M = self.assemble()
        free = self._free_dofs()
        Kr = K[np.ix_(free, free)]
        Mr = M[np.ix_(free, free)]
        am, ak = self._damping(rayleigh)
        Cr = am * Mr + ak * Kr

        dt = t_final / n_steps
        nf = len(free)

        # initial conditions
        u0 = np.zeros(self.n_dof)
        if preload_tip or preload_distributed:
            u0 = self._static_full(preload_tip, preload_distributed)
        if initial_displacement is not None:
            if preload_tip or preload_distributed:
                raise ValueError("choose preload or initial_displacement, not both")
            u0 = self._vector(initial_displacement, "initial_displacement")
        v0 = np.zeros(self.n_dof) if initial_velocity is None else self._vector(initial_velocity, "initial_velocity")
        if np.any(u0[self.constrained_dofs] != 0) or np.any(v0[self.constrained_dofs] != 0):
            raise ValueError("initial conditions must vanish at constrained DOFs")
        u = u0[free].copy()
        v = v0[free].copy()

        def fvec(t):
            return self._vector(force(t), "force(t)")[free] if force is not None else np.zeros(nf)

        a = np.linalg.solve(Mr, fvec(0.0) - Cr @ v - Kr @ u)

        # Newmark constants
        a0 = 1.0 / (beta * dt ** 2)
        a1 = gamma / (beta * dt)
        a2 = 1.0 / (beta * dt)
        a3 = 1.0 / (2 * beta) - 1.0
        a4 = gamma / beta - 1.0
        a5 = dt * (gamma / (2 * beta) - 1.0)

        Keff = Kr + a0 * Mr + a1 * Cr
        # factorise once (constant dt, linear system)
        from scipy.linalg import lu_factor, lu_solve
        lu = lu_factor(Keff)

        n_out = n_steps + 1
        u_hist = np.zeros((n_out, self.n_dof))
        v_hist = np.zeros_like(u_hist)
        a_hist = np.zeros_like(u_hist)
        t = np.linspace(0.0, t_final, n_out)

        u_hist[0, free], v_hist[0, free], a_hist[0, free] = u, v, a

        for n in range(1, n_out):
            tn = t[n]
            rhs = (fvec(tn)
                   + Mr @ (a0 * u + a2 * v + a3 * a)
                   + Cr @ (a1 * u + a4 * v + a5 * a))
            u_new = lu_solve(lu, rhs)
            a_new = a0 * (u_new - u) - a2 * v - a3 * a
            v_new = v + dt * ((1 - gamma) * a + gamma * a_new)
            u, v, a = u_new, v_new, a_new

            u_hist[n, free], v_hist[n, free], a_hist[n, free] = u, v, a

        kinetic = 0.5 * np.einsum("ij,jk,ik->i", v_hist[:, free], Mr, v_hist[:, free])
        strain = 0.5 * np.einsum("ij,jk,ik->i", u_hist[:, free], Kr, u_hist[:, free])
        return TransientResult(t, u_hist, v_hist, a_hist, kinetic, strain)

    def transient_response(self, t_final, n_steps, preload_tip=0.0,
                           preload_distributed=0.0, force=None,
                           rayleigh=(0.0, 0.0), beta=0.25, gamma=0.5):
        """Compatibility API returning time, right-end w, and all nodal w."""
        result = self.transient_result(t_final, n_steps, preload_tip,
                                       preload_distributed, force, rayleigh, beta, gamma)
        return result.time, result.deflection[:, -1], result.deflection

    def _uniform_reference(self):
        if self.EI_profile is not None or self.rhoA_profile is not None:
            raise ValueError("analytical references require a uniform beam without profiles")

    # --- analytical references (for validation/plots) ---
    def analytical_tip_deflection(self, tip_load=0.0, distributed_load=0.0):
        self._uniform_reference()
        if self.boundary != "cantilever":
            raise ValueError("analytical_tip_deflection applies only to cantilevers")
        EI = self.E * self.I
        L = self.length
        w = 0.0
        if tip_load:
            w += tip_load * L ** 3 / (3 * EI)
        if distributed_load:
            w += distributed_load * L ** 4 / (8 * EI)
        return w

    def analytical_frequencies(self, n_modes=4):
        self._uniform_reference()
        _count(n_modes, "n_modes")
        EI = self.E * self.I
        rhoA = self.rho * self.A
        L = self.length
        if self.boundary == "simply-supported":
            bl = np.arange(1, n_modes + 1) * np.pi
        elif self.boundary == "clamped-clamped":
            roots = np.array([4.7300407449, 7.8532046241, 10.9956078380,
                              14.1371654913, 17.2787596574])
            if n_modes > len(roots):
                raise ValueError("clamped-clamped analytical reference supports at most five modes")
            bl = roots[:n_modes]
        else:
            if n_modes > len(CANTILEVER_BETA_L):
                raise ValueError("cantilever analytical reference supports at most five modes")
            bl = CANTILEVER_BETA_L[:n_modes]
        omega = bl ** 2 * np.sqrt(EI / (rhoA * L ** 4))
        return omega / (2 * np.pi)


def steel_rectangular(length=1.0, width=0.05, height=0.01, n_elements=20, boundary="cantilever"):
    """Convenience constructor: a steel beam with a rectangular cross-section."""
    A, I = rectangular_section(width, height)
    return Beam(length=length, E=200e9, I=I, rho=7850.0, A=A, n_elements=n_elements, boundary=boundary)


# Common engineering materials: (Young's modulus [Pa], density [kg/m^3])
MATERIALS = {
    "steel":    (200e9, 7850.0),
    "aluminum": (69e9,  2710.0),
    "copper":   (117e9, 8960.0),
    "titanium": (116e9, 4506.0),
}


def beam_for_material(material, length=1.0, width=0.05, height=0.01, n_elements=20, boundary="cantilever"):
    """Build a rectangular beam for a named material from ``MATERIALS``."""
    if material not in MATERIALS:
        raise ValueError(f"material must be one of {', '.join(MATERIALS)}")
    E, rho = MATERIALS[material]
    A, I = rectangular_section(width, height)
    return Beam(length=length, E=E, I=I, rho=rho, A=A, n_elements=n_elements, boundary=boundary)


def rectangular_section(width, height):
    """Return (A, I) for a solid rectangular cross-section."""
    _positive(width, "width")
    _positive(height, "height")
    return width * height, width * height ** 3 / 12.0


def rayleigh_coefficients(frequencies, damping_ratios):
    """Fit alpha_m [1/s], alpha_k [s] to two modal damping ratios.

    zeta_i = alpha_m/(2 omega_i) + alpha_k omega_i/2.
    Reject a fit requiring negative coefficients (not positive damping globally).
    """
    f = np.asarray(frequencies, dtype=float)
    zeta = np.asarray(damping_ratios, dtype=float)
    if f.shape != (2,) or not np.all(np.isfinite(f)) or np.any(f <= 0) or f[0] == f[1]:
        raise ValueError("provide two distinct finite positive frequencies")
    if zeta.shape != (2,) or not np.all(np.isfinite(zeta)) or np.any(zeta < 0):
        raise ValueError("provide two finite nonnegative damping ratios")
    omega = 2 * np.pi * f
    coefficients = np.linalg.solve(np.column_stack((1 / (2 * omega), omega / 2)), zeta)
    if np.any(coefficients < -1e-12):
        raise ValueError("damping targets require a negative Rayleigh coefficient")
    return tuple(np.maximum(coefficients, 0))
