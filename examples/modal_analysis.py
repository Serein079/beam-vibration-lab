"""Natural frequencies and mode shapes of a steel cantilever beam, compared
against the analytical Euler-Bernoulli eigenvalues."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import fem_beam as fb  # noqa: E402


def main():
    beam = fb.steel_rectangular(length=1.0, width=0.05, height=0.01, n_elements=60)
    n_modes = 4
    freqs, shapes = beam.natural_frequencies(n_modes)
    fa = beam.analytical_frequencies(n_modes)

    print("Cantilever natural frequencies (steel, L=1 m, 50x10 mm section)")
    print("-" * 60)
    print(f"  {'mode':>4} {'FEM (Hz)':>12} {'analytical (Hz)':>16} {'rel. err':>10}")
    for k in range(n_modes):
        print(f"  {k + 1:>4} {freqs[k]:>12.3f} {fa[k]:>16.3f} "
              f"{abs(freqs[k] - fa[k]) / fa[k]:>10.2e}")

    x = beam.node_x
    fig, ax = plt.subplots(figsize=(8, 5))
    for k in range(n_modes):
        ax.plot(x, shapes[k], label=f"mode {k + 1}: {freqs[k]:.1f} Hz")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("x along beam (m)")
    ax.set_ylabel("normalised mode shape")
    ax.set_title("Cantilever beam mode shapes (free vibration)")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "..", "results", "mode_shapes.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"  saved -> {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
