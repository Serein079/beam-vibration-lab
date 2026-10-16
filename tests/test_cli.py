import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from fem_beam_cli import run_analysis


def case(kind):
    return {"beam": {"length": 1, "width": 0.05, "height": 0.01, "E": 200e9,
                     "rho": 7850, "n_elements": 8}, "analysis": {"type": kind}}


@pytest.mark.parametrize("kind", ["static", "modal", "transient", "harmonic"])
def test_json_exports_and_replay(tmp_path, kind):
    config = case(kind)
    if kind == "static":
        config["analysis"]["nodal_forces"] = {"8": -1}
    if kind == "transient":
        config["analysis"].update(t_final=0.02, n_steps=100, tip_load=-1)
    if kind == "harmonic":
        config["analysis"].update(frequencies=[0, 5, 10], tip_load=-1, rayleigh=[2, 0])
    result = run_analysis(config, tmp_path)
    assert result["analysis"] == kind
    assert json.loads((tmp_path / "input.json").read_text()) == config
    assert (tmp_path / "data.csv").read_text().splitlines()[0]
    with np.load(tmp_path / "solution.npz") as data:
        assert "x" in data
        assert len(data["x"]) == 9
    run_analysis(json.loads((tmp_path / "input.json").read_text()), tmp_path / "replay")
    assert json.loads((tmp_path / "replay" / "summary.json").read_text()) == result


def test_cli_subprocess_and_typo_rejection(tmp_path):
    config = case("static")
    config["analysis"]["tip_load"] = -1
    file = tmp_path / "case.json"
    file.write_text(json.dumps(config))
    env_path = str(Path(__file__).resolve().parents[1] / "src")
    import os
    env = dict(os.environ, PYTHONPATH=env_path)
    result = subprocess.run([sys.executable, "-m", "fem_beam_cli", "--config", str(file),
                             "--output", str(tmp_path / "out")], env=env, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout)["right_end_deflection_m"] == pytest.approx(-0.0004, rel=1e-8)
    config["analysis"]["typo"] = 10
    file.write_text(json.dumps(config))
    result = subprocess.run([sys.executable, "-m", "fem_beam_cli", "--config", str(file),
                             "--output", str(tmp_path / "bad")], env=env, capture_output=True, text=True)
    assert result.returncode == 2
    assert "unknown analysis keys" in result.stderr
    assert not (tmp_path / "bad").exists()
