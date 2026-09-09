"""Run clean-room n=15 and n=21 brute-force reduction regressions."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def parse_line(line: str) -> dict:
    fields = line.split()
    if fields[0] != "ORDER":
        raise ValueError(line)
    result = {"order": int(fields[1])}
    for field in fields[2:]:
        key, value = field.split("=", 1)
        if value in {"0", "1"} and key not in {"maximum", "root", "modulus", "brute_words", "normalized_words", "profile_groups", "fold3_signatures", "fold5_signatures", "fold7_signatures", "brute_maximizers", "reduced_maximizers", "brute_affine_classes", "reduced_affine_classes", "brute_affine_hash", "reduced_affine_hash"}:
            result[key] = value == "1"
        elif value.lstrip("-").isdigit():
            result[key] = int(value)
        else:
            result[key] = value
    return result


def run(output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    source = Path(__file__).with_suffix(".cpp")
    exe = output / "small_order_validation.exe"
    subprocess.run(["clang++", "-O3", "-std=c++20", "-Wall", "-Wextra", str(source), "-o", str(exe)], check=True)
    reports = []
    for order, factor in ((15, 5), (21, 7)):
        completed = subprocess.run([str(exe), str(order), str(factor)], capture_output=True, text=True, check=True)
        reports.append(parse_line(completed.stdout.strip()))
    report = {
        "status": "PASS" if all(row["status"] == "PASS" for row in reports) else "FAIL",
        "implementation": "audit/completeness/small_order_validation.cpp",
        "production_executables_called": False,
        "production_verifier_called": False,
        "orders": {str(row["order"]): row for row in reports},
        "interpretation": "This is an end-to-end regression of complement normalization, factor folds, profile grouping, and exact fiber reconstruction against brute force; it is not an independent order-55 fiber-union proof.",
    }
    (output / "small_order_validation_report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    for row in reports:
        print("ORDER", row["order"], "status=", row["status"], "maximum=", row["maximum"], "profile_groups=", row["profile_groups"], "coverage=", row["coverage"])
    print("SMALL-ORDER END-TO-END:", report["status"])
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("certificates/order55_global_v2"))
    args = parser.parse_args()
    raise SystemExit(run(args.output.resolve()))
