"""Non-uniform (cut-out) beam: a cantilever with periodic lightening holes that
locally reduce the bending stiffness, modelled with a per-element EI profile.

This shows how removing material (e.g. bolt or lightening holes) increases tip
deflection and lowers the natural frequencies, compared against the equivalent
solid beam."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import fem_beam as fb  # noqa: E402


def main():
    length, width, height = 1.0, 0.05, 0.02
    n = 60
    A, I = fb.rectangular_section(width, height)
    EI_solid = 200e9 * I
    rhoA_solid = 7850.0 * A

    # Four weakened regions: 60% local EI reduction, 40% local rhoA reduction.
    # A surrogate profile, not an explicitly meshed hole geometry.
    EI_profile = np.full(n, EI_solid)
    rhoA_profile = np.full(n, rhoA_solid)
    xc = np.linspace(0, n, 6)[1:-1].astype(int)
    for c in xc:
        for e in (c - 1, c, c + 1):
            if 0 <= e < n:
                EI_profile[e] *= 0.4
                rhoA_profile[e] *= 0.6

    solid = fb.Beam(length, 200e9, I, 7850.0, A, n_elements=n)
    holed = fb.Beam(length, 200e9, I, 7850.0, A, n_elements=n,
                    EI_profile=EI_profile, rhoA_profile=rhoA_profile)

    _, w_solid = solid.static_deflection(tip_load=-200.0)
    _, w_holed = holed.static_deflection(tip_load=-200.0)
    f_solid = solid.natural_frequencies(3)[0]
    f_holed = holed.natural_frequencies(3)[0]

    print("Solid vs cut-out cantilever (L=1 m, 200 N tip load)")
    print("-" * 56)
    print(f"  tip deflection : solid {w_solid[-1]*1e3:7.3f} mm | "
          f"cut-out {w_holed[-1]*1e3:7.3f} mm")
    print(f"  f1 (Hz)        : solid {f_solid[0]:7.2f}    | cut-out {f_holed[0]:7.2f}")
    print(f"  f2 (Hz)        : solid {f_solid[1]:7.2f}    | cut-out {f_holed[1]:7.2f}")
    print(f"  total mass ratio: {holed.total_mass / solid.total_mass:.3f} (8% reduction)")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    ax1.plot(solid.node_x, w_solid * 1e3, label="solid")
    ax1.plot(holed.node_x, w_holed * 1e3, label="with cut-outs")
    ax1.set_xlabel("x (m)"); ax1.set_ylabel("deflection (mm)")
    ax1.set_title("Static deflection under tip load")
    ax1.grid(True, alpha=0.3); ax1.legend()
    ax2.step(np.linspace(0, length, n), EI_profile / EI_solid, where="mid", color="#228833")
    ax2.set_xlabel("x (m)"); ax2.set_ylabel("EI / EI$_{solid}$")
    ax2.set_title("Bending-stiffness profile (cut-outs)")
    ax2.set_ylim(0, 1.1); ax2.grid(True, alpha=0.3)
    fig.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "..", "results", "stepped_beam_cutouts.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"  saved -> {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
