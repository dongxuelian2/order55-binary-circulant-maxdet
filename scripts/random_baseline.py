"""Run the deliberately small random baseline smoke experiment."""

from __future__ import annotations

import argparse
from datetime import datetime, timezone
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

from maxdet import run_random_baseline


def _git_commit(repo_root: Path) -> str:
    try:
        completed = subprocess.run(
            ["git", "-c", f"safe.directory={repo_root}", "rev-parse", "HEAD"],
            cwd=repo_root,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return completed.stdout.strip()


def _package_versions() -> dict[str, str]:
    versions: dict[str, str] = {}
    for name in ("maximal-determinant", "numpy", "scipy", "sympy", "pytest", "numba"):
        try:
            versions[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            versions[name] = "not installed"
    return versions


def _record(result, repo_root: Path) -> dict[str, object]:
    record = result.as_dict()
    record.update(
        {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "git_commit": _git_commit(repo_root),
            "python_version": platform.python_version(),
            "packages": _package_versions(),
        }
    )
    return record


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n", type=int, required=True, help="matrix order")
    parser.add_argument("--samples", type=int, required=True, help="number of random candidates")
    parser.add_argument("--seed", type=int, required=True, help="NumPy RNG seed")
    parser.add_argument(
        "--output",
        type=Path,
        help="optional JSON output path; metadata is included when writing a file",
    )
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    result = run_random_baseline(args.n, args.samples, args.seed)
    record = _record(result, repo_root)

    print(f"n: {record['n']}")
    print(f"samples: {record['samples']}")
    print(f"seed: {record['seed']}")
    print(f"best fast score: {record['best_fast_score']:.17g}")
    print(f"exact determinant: {record['exact_determinant']}")
    print("matrix:")
    for row in record["matrix"]:
        print("  " + " ".join(f"{value:2d}" for value in row))

    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(
            json.dumps(record, indent=2) + "\n",
            encoding="utf-8",
        )
        print(f"wrote: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
