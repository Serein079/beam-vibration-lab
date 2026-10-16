"""Transient dynamics: a cantilever is deflected statically by a tip load, then
released and allowed to vibrate freely. The response is integrated in time with
the Newmark-beta method, with and without damping.

The free-vibration (ring-down) frequency should match the fundamental natural
frequency from the modal analysis -- a strong cross-check of the two solvers."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import fem_beam as fb  # noqa: E402


def main():
    beam = fb.steel_rectangular(length=1.0, width=0.05, height=0.01, n_elements=40)
    f1 = beam.natural_frequencies(1)[0][0]

    t_final, n_steps = 0.6, 6000
    t, tip_u, _ = beam.transient_response(t_final, n_steps, preload_tip=-100.0,
                                          rayleigh=(0.0, 0.0))
    t_d, tip_d, _ = beam.transient_response(t_final, n_steps, preload_tip=-100.0,
                                            rayleigh=(2.0, 2e-5))

    # ring-down frequency from zero crossings of the undamped response
    s = np.sign(tip_u - tip_u.mean())
    crossings = t[np.where(np.diff(s) != 0)[0]]
    f_ring = 1.0 / (2 * np.mean(np.diff(crossings)))

    print("Transient ring-down of a released cantilever (steel, L=1 m)")
    print("-" * 60)
    print(f"  modal fundamental frequency : {f1:.3f} Hz")
    print(f"  measured ring-down frequency: {f_ring:.3f} Hz  "
          f"(rel. diff {abs(f_ring - f1) / f1:.1e})")
    print(f"  initial tip deflection      : {tip_u[0] * 1e3:.2f} mm")

    fig, ax = plt.subplots(figsize=(9, 4.6))
    ax.plot(t, tip_u * 1e3, lw=1.0, label="undamped")
    ax.plot(t_d, tip_d * 1e3, lw=1.4, label="Rayleigh damping")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("time (s)")
    ax.set_ylabel("tip deflection (mm)")
    ax.set_title(f"Cantilever free-vibration ring-down  "
                 f"($f_1$ = {f1:.2f} Hz, Newmark-$\\beta$)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "..", "results", "transient_ringdown.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"  saved -> {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
