"""Create a deterministic source ZIP without environments or execution caches."""
import argparse
import hashlib
from pathlib import Path
import zipfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("dist"))
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    args.output.mkdir(parents=True, exist_ok=True)
    destination = args.output / "beam-vibration-lab-0.2.0-source.zip"
    roots = {"src", "examples", "tests", "configs", "docs", "fenics", "scripts", "results", ".github"}
    files = [p for p in root.rglob("*") if p.is_file()
             and (len(p.relative_to(root).parts) == 1 or p.relative_to(root).parts[0] in roots)
             and "__pycache__" not in p.parts
             and not any(part.endswith(".egg-info") for part in p.parts)
             and p.suffix not in (".pyc", ".pyo")]
    with zipfile.ZipFile(destination, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for path in sorted(files):
            info = zipfile.ZipInfo("beam-vibration-lab/" + path.relative_to(root).as_posix(),
                                   date_time=(2026, 10, 9, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes(), compresslevel=9)
    digest = hashlib.sha256(destination.read_bytes()).hexdigest()
    destination.with_suffix(".zip.sha256").write_text(f"{digest}  {destination.name}\n", encoding="ascii")
    print(f"{destination}: {len(files)} files, SHA256 {digest}")


if __name__ == "__main__":
    main()
