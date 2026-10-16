"""Upstream regression tests validating the FEM beam solver against
closed-form analytical solutions.

Run:  python -m pytest    (or: python tests/test_fem_beam.py)
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import fem_beam as fb  # noqa: E402


def _beam(n):
    return fb.steel_rectangular(length=1.0, width=0.05, height=0.01, n_elements=n)


def test_tip_point_load_matches_PL3_over_3EI():
    beam = _beam(20)
    _, w = beam.static_deflection(tip_load=-100.0)
    exact = -beam.analytical_tip_deflection(tip_load=100.0)
    # Cubic elements are exact for a tip point load.
    assert np.isclose(w[-1], exact, rtol=1e-6)


def test_distributed_load_matches_qL4_over_8EI():
    beam = _beam(40)
    _, w = beam.static_deflection(distributed_load=-200.0)
    exact = -beam.analytical_tip_deflection(distributed_load=200.0)
    assert np.isclose(w[-1], exact, rtol=1e-4)


def test_natural_frequencies_match_euler_bernoulli():
    beam = _beam(60)
    f = beam.natural_frequencies(4)[0]
    fa = beam.analytical_frequencies(4)
    assert np.allclose(f, fa, rtol=2e-3)


def test_clamped_end_has_zero_deflection_and_slope():
    beam = _beam(10)
    _, w = beam.static_deflection(tip_load=-50.0)
    assert np.isclose(w[0], 0.0, atol=1e-12)
    assert beam.static_result(tip_load=-50.0).rotation[0] == 0.0


def test_convergence_improves_with_refinement():
    coarse = _beam(2).natural_frequencies(1)[0][0]
    fine = _beam(64).natural_frequencies(1)[0][0]
    exact = _beam(2).analytical_frequencies(1)[0]
    assert abs(fine - exact) < abs(coarse - exact)


def test_transient_ringdown_frequency_matches_modal():
    beam = _beam(40)
    f1 = beam.natural_frequencies(1)[0][0]
    t, tip, _ = beam.transient_response(t_final=0.5, n_steps=4000,
                                        preload_tip=-100.0)
    s = np.sign(tip - tip.mean())
    crossings = t[np.where(np.diff(s) != 0)[0]]
    f_ring = 1.0 / (2 * np.mean(np.diff(crossings)))
    assert abs(f_ring - f1) / f1 < 0.02


def test_transient_initial_amplitude_equals_static():
    beam = _beam(30)
    _, w = beam.static_deflection(tip_load=-100.0)
    _, tip, _ = beam.transient_response(t_final=0.1, n_steps=500,
                                        preload_tip=-100.0)
    assert np.isclose(tip[0], w[-1], rtol=1e-9)


def test_damping_decays_amplitude():
    beam = _beam(30)
    _, tip, _ = beam.transient_response(t_final=0.6, n_steps=4000,
                                        preload_tip=-100.0, rayleigh=(3.0, 3e-5))
    early = np.abs(tip[:1000]).max()
    late = np.abs(tip[-1000:]).max()
    assert late < 0.6 * early


def test_uniform_profile_equals_scalar_section():
    plain = _beam(10)
    K, M = plain.assemble()
    n = 10
    prof = fb.Beam(plain.length, plain.E, plain.I, plain.rho, plain.A,
                   n_elements=n,
                   EI_profile=np.full(n, plain.E * plain.I),
                   rhoA_profile=np.full(n, plain.rho * plain.A))
    Kp, Mp = prof.assemble()
    assert np.allclose(K, Kp) and np.allclose(M, Mp)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"PASS {name}")
    print("All FEM regression tests passed.")
