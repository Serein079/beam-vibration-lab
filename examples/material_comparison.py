"""Material study: how the choice of material changes the static deflection and
natural frequencies of an identical beam (steel / aluminum / copper / titanium).

A nice result falls out: steel and aluminum have almost the same fundamental
frequency, because frequency scales with the *specific* stiffness sqrt(E/rho),
which is nearly equal for the two metals."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import fem_beam as fb  # noqa: E402

MATERIALS = ["steel", "aluminum", "copper", "titanium"]


def main():
    load = -50.0  # N tip load
    rows = []
    for m in MATERIALS:
        beam = fb.beam_for_material(m, length=0.5, width=0.02, height=0.005, n_elements=40)
        _, w = beam.static_deflection(tip_load=load)
        f = beam.natural_frequencies(3)[0]
        rows.append((m, w[-1] * 1e3, f))

    print("Material comparison (L=0.5 m, 20x5 mm section, 50 N tip load)")
    print("-" * 64)
    print(f"  {'material':>9} {'tip defl (mm)':>14} {'f1 (Hz)':>9} {'f2 (Hz)':>9} {'f3 (Hz)':>9}")
    for m, tip, f in rows:
        print(f"  {m:>9} {tip:>14.3f} {f[0]:>9.1f} {f[1]:>9.1f} {f[2]:>9.1f}")

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.5))
    names = [r[0] for r in rows]
    tips = [abs(r[1]) for r in rows]
    f1s = [r[2][0] for r in rows]
    ax1.bar(names, tips, color="#4477aa")
    ax1.set_ylabel("|tip deflection| (mm)")
    ax1.set_title("Static deflection under 50 N tip load")
    ax1.grid(True, axis="y", alpha=0.3)
    ax2.bar(names, f1s, color="#cc6677")
    ax2.set_ylabel("fundamental frequency $f_1$ (Hz)")
    ax2.set_title("First natural frequency")
    ax2.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "..", "results", "material_comparison.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"  saved -> {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
