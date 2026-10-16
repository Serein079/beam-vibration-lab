"""Compare supports, plot a frequency response, and check Newmark energy."""
import json
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from fem_beam import steel_rectangular, rayleigh_coefficients


def main():
    out = Path(__file__).resolve().parents[1] / "results"
    out.mkdir(exist_ok=True)
    supports = ["cantilever", "simply-supported", "clamped-clamped"]
    summary = {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    for support in supports:
        b = steel_rectangular(n_elements=30, boundary=support)
        r = b.static_result(distributed_load=-20)
        modal = b.modal_result(3)
        summary[support] = {"max_deflection_mm": float(np.max(np.abs(r.deflection)) * 1000),
                            "frequencies_Hz": modal.frequencies.tolist()}
        axes[0].plot(r.x, r.deflection * 1000, label=support)
        axes[1].plot([1, 2, 3], modal.frequencies, "o-", label=support)
    axes[0].set(xlabel="x (m)", ylabel="deflection (mm)", title="Supports: uniform load -20 N/m")
    axes[1].set(xlabel="mode", ylabel="frequency (Hz)", title="Support-dependent natural frequencies", xticks=[1, 2, 3])
    for ax in axes:
        ax.grid(alpha=.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(out / "support_comparison.png", dpi=150)
    plt.close(fig)

    b = steel_rectangular(n_elements=30)
    f = b.modal_result(2).frequencies
    damping = rayleigh_coefficients(f, [0.02, 0.02])
    frequencies = np.linspace(0, 1.8 * f[1], 800)
    response = b.harmonic_response(frequencies, b.load_vector(tip_load=-1), damping)
    fig, axes = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    axes[0].semilogy(frequencies, np.abs(response[:, -2]) * 1000)
    axes[0].set(ylabel="tip amplitude (mm)", title="Tip harmonic force 1 N, 2% damping at modes 1 and 2")
    axes[1].plot(frequencies, np.angle(response[:, -2], deg=True))
    axes[1].set(xlabel="excitation frequency (Hz)", ylabel="phase (deg)")
    for ax in axes:
        for frequency in f:
            ax.axvline(frequency, color="gray", ls="--", alpha=.5)
        ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(out / "harmonic_response.png", dpi=150)
    plt.close(fig)
    np.savetxt(out / "harmonic_response.csv", np.column_stack((frequencies, np.abs(response[:, -2]), np.angle(response[:, -2]))),
               delimiter=",", header="frequency_Hz,amplitude_m,phase_rad", comments="")

    undamped = b.transient_result(.5, 2000, preload_tip=-1)
    damped = b.transient_result(.5, 2000, preload_tip=-1, rayleigh=damping)
    summary["dynamics"] = {"rayleigh_coefficients": list(damping),
                           "undamped_relative_energy_drift": float(np.max(np.abs(undamped.total_energy / undamped.total_energy[0] - 1))),
                           "damped_final_energy_ratio": float(damped.total_energy[-1] / damped.total_energy[0])}
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].plot(undamped.time, undamped.deflection[:, -1] * 1000, label="undamped")
    axes[0].plot(damped.time, damped.deflection[:, -1] * 1000, label="damped")
    axes[0].set(xlabel="time (s)", ylabel="tip deflection (mm)", title="Released static preload -1 N")
    axes[1].plot(undamped.time, undamped.total_energy / undamped.total_energy[0], label="undamped")
    axes[1].plot(damped.time, damped.total_energy / damped.total_energy[0], label="damped")
    axes[1].set(xlabel="time (s)", ylabel="E(t) / E(0)", title="Mechanical energy")
    for ax in axes:
        ax.grid(alpha=.3)
        ax.legend()
    fig.tight_layout()
    fig.savefig(out / "energy_check.png", dpi=150)
    plt.close(fig)
    (out / "extended_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
