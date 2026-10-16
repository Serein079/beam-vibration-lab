"""JSON-configured analyses; all quantities use SI units."""
import argparse
import json
from pathlib import Path

import numpy as np

from fem_beam import Beam, rectangular_section, __version__


def _keys(mapping, allowed, name):
    if not isinstance(mapping, dict):
        raise ValueError(f"{name} must be a JSON object")
    unknown = set(mapping) - set(allowed)
    if unknown:
        raise ValueError(f"unknown {name} keys: {', '.join(sorted(unknown))}")


def _nodes(mapping):
    if not isinstance(mapping, dict):
        raise ValueError("nodal loads must be objects mapping node indices to values")
    result = {}
    for node, value in mapping.items():
        index = int(node)
        if str(index) != str(node):
            raise ValueError("node keys must be canonical integer strings")
        result[index] = value
    return result


def run_analysis(config, output):
    """Validate a case and export CSV, NPZ and JSON (without plotting dependency)."""
    _keys(config, ("beam", "analysis"), "config")
    data = dict(config["beam"])
    _keys(data, ("length", "E", "rho", "A", "I", "width", "height", "n_elements",
                 "boundary", "EI_profile", "rhoA_profile"), "beam")
    if "width" in data or "height" in data:
        if "A" in data or "I" in data:
            raise ValueError("choose width/height or A/I")
        data["A"], data["I"] = rectangular_section(data.pop("width"), data.pop("height"))
    beam = Beam(**data)
    analysis = config["analysis"]
    _keys(analysis, ("type", "tip_load", "distributed_load", "nodal_forces", "nodal_moments",
                     "n_modes", "frequencies", "rayleigh", "t_final", "n_steps"), "analysis")
    kind = analysis["type"]
    load_keys = {"tip_load", "distributed_load", "nodal_forces", "nodal_moments"}
    allowed = {
        "static": load_keys,
        "modal": {"n_modes"},
        "harmonic": load_keys | {"frequencies", "rayleigh"},
        "transient": load_keys | {"t_final", "n_steps", "rayleigh"},
    }
    if kind not in allowed:
        raise ValueError("analysis type must be static, modal, harmonic or transient")
    if set(analysis) - allowed[kind] - {"type"}:
        raise ValueError(f"unsupported options for {kind}")
    load = beam.load_vector(analysis.get("tip_load", 0), analysis.get("distributed_load", 0),
                            _nodes(analysis.get("nodal_forces", {})),
                            _nodes(analysis.get("nodal_moments", {})))
    summary = {"version": __version__, "analysis": kind, "boundary": beam.boundary,
               "n_elements": beam.n_elements, "mass_kg": beam.total_mass,
               "units": "SI: m, s, N, Pa, kg; frequency Hz; energy J"}
    if kind == "static":
        r = beam.static_result(load=load)
        table = np.column_stack((r.x, r.deflection, r.rotation, r.reactions[0::2], r.reactions[1::2]))
        headers = ["x_m", "w_m", "theta_rad", "reaction_N", "reaction_Nm"]
        arrays = {"x": r.x, "displacement": r.displacement, "reactions": r.reactions, "load": r.load}
        summary["max_abs_deflection_m"] = float(np.max(np.abs(r.deflection)))
        summary["right_end_deflection_m"] = float(r.deflection[-1])
    elif kind == "modal":
        r = beam.modal_result(analysis.get("n_modes", 4))
        table = np.column_stack((np.arange(1, len(r.frequencies) + 1), r.frequencies))
        headers = ["mode", "frequency_Hz"]
        arrays = {"x": beam.node_x, "frequencies": r.frequencies, "modes": r.modes}
        summary["frequencies_Hz"] = r.frequencies.tolist()
    elif kind == "harmonic":
        frequencies = np.asarray(analysis["frequencies"], dtype=float)
        r = beam.harmonic_response(frequencies, load, analysis.get("rayleigh", (0, 0)))
        # Monitor a free translational DOF, since right-end w can be constrained.
        w_dofs = beam._free_dofs()[beam._free_dofs() % 2 == 0]
        if len(w_dofs) == 0:
            raise ValueError("harmonic CSV needs at least one free translational DOF")
        monitor = beam.n_dof - 2 if beam.n_dof - 2 in w_dofs else int(w_dofs[len(w_dofs) // 2])
        table = np.column_stack((frequencies, np.abs(r[:, monitor]), np.angle(r[:, monitor]),
                                 r[:, monitor].real, r[:, monitor].imag))
        headers = ["frequency_Hz", "amplitude_m", "phase_rad", "real_m", "imag_m"]
        arrays = {"x": beam.node_x, "frequencies": frequencies, "response": r}
        summary["monitor_node"] = monitor // 2
    else:
        initial = beam.static_result(load=load).displacement
        r = beam.transient_result(analysis["t_final"], analysis["n_steps"],
                                  initial_displacement=initial,
                                  rayleigh=analysis.get("rayleigh", (0, 0)))
        table = np.column_stack((r.time, r.deflection[:, -1],
                                 np.max(np.abs(r.deflection), axis=1),
                                 r.kinetic_energy, r.strain_energy, r.total_energy))
        headers = ["time_s", "right_end_w_m", "max_abs_w_m", "kinetic_J", "strain_J", "total_J"]
        arrays = {"x": beam.node_x, "time": r.time, "displacement": r.displacement,
                  "velocity": r.velocity, "acceleration": r.acceleration,
                  "kinetic_energy": r.kinetic_energy, "strain_energy": r.strain_energy}
        summary["transient_kind"] = "static preload released at t=0 (free ring-down)"
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    np.savetxt(output / "data.csv", table, delimiter=",", header=",".join(headers), comments="")
    np.savez_compressed(output / "solution.npz", **arrays)
    (output / "summary.json").write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (output / "input.json").write_text(json.dumps(config, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=__version__)
    parser.add_argument("--config", required=True, type=Path, help="JSON case file")
    parser.add_argument("--output", required=True, type=Path, help="result directory")
    args = parser.parse_args(argv)
    try:
        config = json.loads(args.config.read_text(encoding="utf-8"))
        summary = run_analysis(config, args.output)
    except (ValueError, TypeError, KeyError, OSError, np.linalg.LinAlgError) as exc:
        parser.exit(2, f"beam-fem: {exc}\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
