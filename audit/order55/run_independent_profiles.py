"""Run and compare the clean-room correlation-profile screen."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CERT = ROOT / "certificates" / "order55_global"
SOURCE = ROOT / "audit" / "order55" / "independent_profile_screen.cpp"
MAXIMUM = "134694094094758395331307111329132"
PARTITIONS = ((25, 2), (26, 3), (26, 5), (26, 6),
              (27, 1), (27, 3), (27, 4), (27, 5))


def read_jsonl(path: Path) -> list[dict]:
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def target_key(row: dict) -> tuple[str, tuple[int, ...]]:
    return (str(row["absolute_profile_product"]), tuple(row["correlations"]))


def run_one(executable: Path, k: int, partition: int, directory: Path) -> dict:
    stem = f"corr_k{k}_part{partition}"
    actual_audit = directory / f"{stem}.audit.json"
    actual_targets = directory / f"{stem}.targets.jsonl"
    completed = subprocess.run(
        [str(executable), str(k), MAXIMUM,
         str(CERT / f"{stem}.txt"), str(CERT),
         str(actual_audit), str(actual_targets)],
        capture_output=True, text=True, check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"{stem} failed with exit {completed.returncode}: "
                           f"{completed.stdout}\n{completed.stderr}")
    actual_audit_row = json.loads(actual_audit.read_text().strip())
    expected_audit = json.loads((CERT / f"{stem}_audit.json").read_text())
    observed = {
        key: actual_audit_row[key]
        for key in ("k", "ordered_profiles", "folded_matching_profiles",
                    "canonical_determinants", "above_screen")
    }
    if observed != {
        key: expected_audit[key]
        for key in ("k", "ordered_profiles", "folded_matching_profiles",
                    "canonical_determinants", "above_screen")
    }:
        raise AssertionError((stem, "counter mismatch", actual_audit_row, expected_audit))
    actual = read_jsonl(actual_targets)
    expected = read_jsonl(CERT / f"{stem}_targets.jsonl")
    if {target_key(row) for row in actual} != {target_key(row) for row in expected}:
        raise AssertionError((stem, "target set mismatch", len(actual), len(expected)))
    return {
        "k": k,
        "partition": partition,
        "ordered_profiles": actual_audit_row["ordered_profiles"],
        "folded_matching_profiles": actual_audit_row["folded_matching_profiles"],
        "canonical_determinants": actual_audit_row["canonical_determinants"],
        "above_screen": actual_audit_row["above_screen"],
        "max_absolute_profile_product": actual_audit_row["max_absolute_profile_product"],
        "targets": len(actual),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    args = parser.parse_args()
    with tempfile.TemporaryDirectory(prefix="order55-profile-audit-") as temporary:
        directory = Path(temporary)
        executable = directory / "independent_profile_screen.exe"
        compile_result = subprocess.run(
            ["clang++", "-O3", "-std=c++20", "-Wall", "-Wextra",
             "-Wconversion", "-Wsign-conversion", str(SOURCE), "-o", str(executable)],
            capture_output=True, text=True, check=False,
        )
        if compile_result.returncode != 0:
            raise RuntimeError(f"compile failed:\n{compile_result.stdout}\n{compile_result.stderr}")
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            futures = [pool.submit(run_one, executable, k, part, directory)
                       for k, part in PARTITIONS]
            rows = [future.result() for future in futures]
    rows.sort(key=lambda row: (row["k"], row["partition"]))
    summary = {
        "partitions": rows,
        "targets": sum(row["targets"] for row in rows),
        "status": "PASS",
    }
    output = ROOT / "audit" / "order55" / "independent_profile_audit.json"
    output.write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, separators=(",", ":")))


if __name__ == "__main__":
    main()
