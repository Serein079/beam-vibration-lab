"""Independent physical checks for the extensions, plus invalid-input checks."""
import numpy as np
import pytest

from fem_beam import Beam, steel_rectangular, rayleigh_coefficients, rectangular_section


def test_simply_supported_uniform_load_and_reactions():
    b = steel_rectangular(n_elements=20, boundary="simply-supported")
    r = b.static_result(distributed_load=-200)
    exact = -5 * 200 * b.length**4 / (384 * b.E * b.I)
    assert r.deflection[10] == pytest.approx(exact, rel=1e-8)
    assert r.reactions[0] == pytest.approx(100, rel=1e-8)
    assert r.reactions[-2] == pytest.approx(100, rel=1e-8)
    assert np.all(r.deflection[[0, -1]] == 0)


def test_fixed_fixed_uniform_load_and_end_moments():
    b = steel_rectangular(n_elements=20, boundary="clamped-clamped")
    r = b.static_result(distributed_load=-200)
    assert r.deflection[10] == pytest.approx(-200 * b.length**4 / (384 * b.E * b.I), rel=1e-8)
    assert r.reactions[0] == pytest.approx(100, rel=1e-8)
    assert r.reactions[-2] == pytest.approx(100, rel=1e-8)
    assert r.reactions[1] == pytest.approx(200 / 12, rel=1e-8)
    assert r.reactions[-1] == pytest.approx(-200 / 12, rel=1e-8)
    assert np.all(r.displacement[b.constrained_dofs] == 0)


def test_midpoint_force_simply_supported():
    b = steel_rectangular(n_elements=20, boundary="simply-supported")
    r = b.static_result(load=b.load_vector(nodal_forces={10: -100}))
    assert r.deflection[10] == pytest.approx(-100 * b.length**3 / (48 * b.E * b.I), rel=1e-8)
    assert r.reactions[0] == pytest.approx(50, rel=1e-8)


def test_cantilever_tip_moment_rotation_and_reaction():
    b = steel_rectangular(n_elements=12)
    r = b.static_result(load=b.load_vector(nodal_moments={12: 10}))
    assert r.deflection[-1] == pytest.approx(10 * b.length**2 / (2 * b.E * b.I), rel=1e-8)
    assert r.rotation[-1] == pytest.approx(10 * b.length / (b.E * b.I), rel=1e-8)
    assert r.reactions[1] == pytest.approx(-10, rel=1e-8)
    assert r.rotation[0] == 0


@pytest.mark.parametrize("boundary", ["cantilever", "simply-supported", "clamped-clamped"])
def test_modes_eigen_residual_mass_orthogonality_and_closed_form(boundary):
    b = steel_rectangular(n_elements=30, boundary=boundary)
    r = b.modal_result(4)
    K, M = b.assemble()
    assert np.allclose(r.frequencies, b.analytical_frequencies(4), rtol=5e-5)
    assert np.allclose(r.modes @ M @ r.modes.T, np.eye(4), atol=2e-10)
    free = b._free_dofs()
    for f, phi in zip(r.frequencies, r.modes):
        residual = (K @ phi - (2 * np.pi * f)**2 * M @ phi)[free]
        assert np.linalg.norm(residual) / np.linalg.norm((K @ phi)[free]) < 2e-6


def test_harmonic_zero_frequency_matches_static_and_complex_equilibrium():
    b = steel_rectangular(n_elements=8)
    load = b.load_vector(tip_load=-1).astype(complex) * (1 + 0.3j)
    frequencies = [0, 12, 30]
    u = b.harmonic_response(frequencies, load, rayleigh=(2, 2e-5))
    assert np.allclose(u[0], b.static_result(tip_load=-1).displacement * (1 + 0.3j))
    K, M = b.assemble()
    C = 2 * M + 2e-5 * K
    free = b._free_dofs()
    for f, phi in zip(frequencies, u):
        w = 2 * np.pi * f
        assert np.allclose(((K - w*w*M + 1j*w*C) @ phi)[free], load[free], atol=1e-8)


def test_harmonic_resonance_first_mode():
    b = steel_rectangular(n_elements=10)
    mode = b.modal_result(1)
    K, M = b.assemble()
    f = mode.frequencies[0]
    # A mass-proportional modal force excites only this mode.
    load = M @ mode.modes[0]
    am = 2.0
    u = b.harmonic_response([f], load, rayleigh=(am, 0))[0]
    expected = mode.modes[0] / (1j * am * 2 * np.pi * f)
    assert np.allclose(u, expected, rtol=2e-6, atol=1e-10)


def test_newmark_conserves_energy_for_undamped_free_vibration():
    b = steel_rectangular(n_elements=8)
    r = b.transient_result(0.3, 900, preload_tip=-1)
    assert np.max(np.abs(r.total_energy / r.total_energy[0] - 1)) < 2e-8
    assert np.all(r.displacement[:, b.constrained_dofs] == 0)
    assert np.all(r.velocity[:, b.constrained_dofs] == 0)


def test_damped_free_energy_is_nonincreasing():
    b = steel_rectangular(n_elements=8)
    r = b.transient_result(0.4, 1000, preload_tip=-1, rayleigh=(3, 3e-5))
    assert np.max(np.diff(r.total_energy)) < 1e-12
    assert r.total_energy[-1] < 0.4 * r.total_energy[0]


def test_modal_initial_velocity_matches_single_oscillator():
    b = steel_rectangular(n_elements=6, boundary="simply-supported")
    modal = b.modal_result(1)
    omega = 2 * np.pi * modal.frequencies[0]
    r = b.transient_result(0.01, 1000, initial_velocity=modal.modes[0])
    exact = np.sin(omega * r.time)[:, None] * modal.modes[0] / omega
    assert np.max(np.abs(r.displacement - exact)) < 1e-7
    assert np.all(r.deflection[:, -1] == 0)


def test_forced_solution_remains_in_static_equilibrium():
    b = steel_rectangular(n_elements=5)
    load = b.load_vector(tip_load=-1)
    initial = b.static_result(load=load).displacement
    r = b.transient_result(0.1, 100, initial_displacement=initial, force=lambda t: load)
    assert np.allclose(r.displacement, initial, atol=1e-11)


def test_rayleigh_fit_and_total_mass():
    f = np.array([10, 100])
    am, ak = rayleigh_coefficients(f, [0.02, 0.03])
    omega = 2 * np.pi * f
    assert np.allclose(am / (2 * omega) + ak * omega / 2, [0.02, 0.03])
    b = steel_rectangular(n_elements=4)
    b.rhoA_profile = np.array([1, 2, 3, 4])
    assert b.total_mass == pytest.approx(2.5)


@pytest.mark.parametrize("kwargs", [
    {"length": 0}, {"E": -1}, {"rho": np.nan}, {"n_elements": 0},
    {"n_elements": 1.5}, {"n_elements": True}, {"boundary": "unknown"},
    {"boundary": "clamped-clamped", "n_elements": 1},
    {"EI_profile": [1, 2]}, {"rhoA_profile": np.zeros(20)},
])
def test_invalid_beams_are_rejected(kwargs):
    values = dict(length=1, E=200e9, I=1e-8, rho=7850, A=1e-4)
    values.update(kwargs)
    with pytest.raises(ValueError):
        Beam(**values)


@pytest.mark.parametrize("kwargs", [
    {"n_steps": 0}, {"n_steps": 2.5}, {"t_final": 0},
    {"rayleigh": [-1, 0]}, {"beta": 0}, {"beta": 0.1},
    {"gamma": 0.4}, {"force": lambda t: np.zeros(2)},
    {"initial_velocity": np.ones(42)},
])
def test_invalid_transient_parameters_are_rejected(kwargs):
    options = dict(t_final=0.1, n_steps=10)
    options.update(kwargs)
    with pytest.raises(ValueError):
        steel_rectangular().transient_result(**options)


def test_reference_mode_and_load_errors():
    b = steel_rectangular(n_elements=2)
    for n in (0, 5, 1.5):
        with pytest.raises(ValueError):
            b.modal_result(n)
    with pytest.raises(ValueError):
        b.analytical_frequencies(6)
    with pytest.raises(ValueError):
        b.load_vector(nodal_forces={-1: 10})
    with pytest.raises(ValueError):
        b.load_vector(nodal_moments={0: np.nan})
    with pytest.raises(ValueError):
        b.harmonic_response([-1], np.zeros(b.n_dof))
    with pytest.raises(ValueError):
        rectangular_section(-1, 1)
