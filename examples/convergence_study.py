"""Mesh-convergence study: how the FEM error in the fundamental natural
frequency decays as the number of elements increases."""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "src"))
import fem_beam as fb  # noqa: E402


def main():
    counts = [1, 2, 4, 8, 16, 32, 64, 128]
    errors = []
    for n in counts:
        beam = fb.steel_rectangular(length=1.0, width=0.05, height=0.01, n_elements=n)
        f1 = beam.natural_frequencies(1)[0][0]
        f1_exact = beam.analytical_frequencies(1)[0]
        errors.append(abs(f1 - f1_exact) / f1_exact)

    print("Convergence of the fundamental frequency")
    print("-" * 44)
    print(f"  {'elements':>9} {'rel. error':>14}")
    for n, e in zip(counts, errors):
        print(f"  {n:>9} {e:>14.3e}")

    fig, ax = plt.subplots(figsize=(7.5, 5))
    ax.loglog(counts, errors, "o-", label="FEM error in $f_1$")
    # reference slope (cubic elements -> ~4th order in frequency)
    ref = np.array(errors[1]) * (np.array(counts[1]) / np.array(counts, dtype=float)) ** 4
    ax.loglog(counts, ref, "k--", alpha=0.6, label="4th-order reference slope")
    ax.set_xlabel("number of elements")
    ax.set_ylabel("relative error in $f_1$")
    ax.set_title("Mesh convergence of the fundamental frequency")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "..", "results", "convergence.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"  saved -> {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
