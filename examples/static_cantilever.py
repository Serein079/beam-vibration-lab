"""Static deflection of a steel cantilever beam under a tip point load and a
uniform distributed load, compared against the analytical solutions."""
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

    x_tip, w_tip = beam.static_deflection(tip_load=-100.0)         # 100 N down
    x_q, w_q = beam.static_deflection(distributed_load=-200.0)     # 200 N/m down

    a_tip = -beam.analytical_tip_deflection(tip_load=100.0)
    a_q = -beam.analytical_tip_deflection(distributed_load=200.0)

    print("Static cantilever deflection (steel, L=1 m, 50x10 mm section)")
    print("-" * 60)
    print(f"  tip point load 100 N : FEM tip = {w_tip[-1] * 1e3:7.4f} mm "
          f"| analytical = {a_tip * 1e3:7.4f} mm")
    print(f"  uniform load 200 N/m : FEM tip = {w_q[-1] * 1e3:7.4f} mm "
          f"| analytical = {a_q * 1e3:7.4f} mm")

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(x_tip, w_tip * 1e3, "o-", ms=3, label="FEM: 100 N tip load")
    ax.plot(x_q, w_q * 1e3, "s-", ms=3, label="FEM: 200 N/m distributed")
    ax.axhline(0, color="k", lw=0.6)
    ax.plot([1.0], [a_tip * 1e3], "r*", ms=14, label="analytical tip (point)")
    ax.plot([1.0], [a_q * 1e3], "m*", ms=14, label="analytical tip (dist.)")
    ax.set_xlabel("x along beam (m)")
    ax.set_ylabel("deflection w (mm)")
    ax.set_title("Cantilever static deflection: FEM vs analytical")
    ax.grid(True, alpha=0.3)
    ax.legend()
    fig.tight_layout()

    out = os.path.join(os.path.dirname(__file__), "..", "results", "static_deflection.png")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    fig.savefig(out, dpi=130)
    print(f"  saved -> {os.path.relpath(out)}")


if __name__ == "__main__":
    main()
