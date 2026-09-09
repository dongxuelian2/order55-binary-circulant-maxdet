"""Run the deterministic one-million-sample forward coverage audit."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path


def run(certificate: Path, output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    targets = json.loads((certificate / "lift_targets.json").read_text(encoding="utf-8"))
    index = output / "target_profile_index_v8.txt"
    rows = sorted({(entry["target"]["k"], tuple(entry["target"]["correlations"])) for entry in targets})
    index.write_text("".join(f"{k} {' '.join(map(str, profile))}\n" for k, profile in rows), encoding="utf-8", newline="\n")
    bound_file = output / "weight_bounds_v8.txt"
    incumbent = 134694094094758395331307111329132
    bounds = []
    for k in range(25):
        numerator = (55 - k) * (k * (55 - k)) ** 27
        denominator = 54 ** 27
        floor_value = numerator // denominator
        bounds.append({"k": k, "upper_floor": str(floor_value), "pruned": floor_value < incumbent})
    bound_file.write_text("".join(f"{row['k']} {row['upper_floor']}\n" for row in bounds), encoding="utf-8", newline="\n")
    source = Path(__file__).with_name("random_forward_coverage_v6.cpp")
    exe = output / "random_forward_coverage_v6.exe"
    subprocess.run(["clang++", "-O3", "-std=c++20", "-Wall", "-Wextra", str(source), "-o", str(exe)], check=True)
    completed = subprocess.run([str(exe), str(index), str(bound_file)], capture_output=True, text=True, check=True)
    fields = dict(token.split("=", 1) for token in completed.stdout.strip().split())
    report = {
        "status": fields["status"], "implementation": "audit/completeness/random_forward_coverage_v6.cpp",
        "production_executables_called": False, "production_verifier_called": False,
        "sample_count": int(fields["samples"]), "rigorously_pruned": int(fields["rigorously_pruned"]),
        "explicitly_covered": int(fields["explicitly_covered"]), "unknown": int(fields["unknown"]),
        "target_signature_count_after_unit_closure": int(fields["target_signature_count"]),
        "exact_bound_rows": int(fields["bound_rows"]), "exact_bounds": bounds,
        "random_generator": "xorshift64* with fixed seed 0x9e3779b97f4a7c15",
        "classification_rule": "low normalized weight <=24 uses an exact integer Ryser bound; active weights outside the independently audited surviving target-profile index are rigorously pruned, and index members are explicitly covered by a lift task",
    }
    (output / "random_forward_coverage_report_v6.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("RANDOM FORWARD COVERAGE")
    print("samples =", report["sample_count"]); print("RIGOROUSLY PRUNED =", report["rigorously_pruned"])
    print("EXPLICITLY COVERED =", report["explicitly_covered"]); print("UNKNOWN =", report["unknown"])
    print("RANDOM FORWARD COVERAGE:", report["status"])
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=Path("certificates/order55_global"))
    parser.add_argument("--output", type=Path, default=Path("certificates/order55_global_v2"))
    args = parser.parse_args()
    raise SystemExit(run(args.certificate.resolve(), args.output.resolve()))
