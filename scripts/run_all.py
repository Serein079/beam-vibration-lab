"""Run validation, all seven examples, and all four JSON cases; fail on errors."""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import platform
import subprocess
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("runs/latest"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "tmp").mkdir(exist_ok=True)
    env = dict(os.environ, PYTHONPATH=str(root / "src"), MPLBACKEND="Agg",
               MPLCONFIGDIR=str(out / "mplconfig"), PYTHONIOENCODING="utf-8",
               TEMP=str(out / "tmp"), TMP=str(out / "tmp"))
    steps = [("pytest", ["-m", "pytest", "-q", f"--junitxml={out / 'pytest.xml'}"])]
    steps += [(p.stem, [str(p)]) for p in sorted((root / "examples").glob("*.py"))]
    steps += [(p.stem, ["-m", "fem_beam_cli", "--config", str(p), "--output", str(out / p.stem)])
              for p in sorted((root / "configs").glob("*.json"))]
    records = []
    for name, command in steps:
        print(f"Running {name}", flush=True)
        started = dt.datetime.now(dt.timezone.utc)
        result = subprocess.run([sys.executable, *command], cwd=root, env=env,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, encoding="utf-8")
        (out / f"{name}.log").write_text(result.stdout, encoding="utf-8")
        records.append({"name": name, "exit_code": result.returncode,
                        "seconds": (dt.datetime.now(dt.timezone.utc) - started).total_seconds()})
    manifest = {"generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                "python": platform.python_version(), "platform": platform.platform(),
                "passed": all(r["exit_code"] == 0 for r in records),
                "fenics": "NOT RUN: requires separate legacy FEniCS environment", "steps": records}
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    lines = ["# Validation run", "", f"Python: {manifest['python']}",
             f"\nGenerated UTC: {manifest['generated_utc']}", "", "| Step | Status | Seconds |",
             "|---|---|---:|"]
    lines += [f"| {r['name']} | {'PASS' if r['exit_code'] == 0 else 'FAIL'} | {r['seconds']:.2f} |" for r in records]
    lines += ["", "Legacy FEniCS: not executed; outside the 1D solver validation scope.", ""]
    (out / "report.md").write_text("\n".join(lines), encoding="utf-8")
    print(f"{'PASS' if manifest['passed'] else 'FAIL'}: {out / 'report.md'}")
    return 0 if manifest["passed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
