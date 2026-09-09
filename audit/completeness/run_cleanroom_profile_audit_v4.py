"""Run and record the independent order-55 correlation-profile audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

M = 134694094094758395331307111329132


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def run(cert: Path, output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    cases = load(cert / "correlation_profiles.json")
    case_file = output / "profile_cases.txt"
    case_file.write_text(
        "".join(f"{i} {row['k']} {' '.join(map(str, row['values']))}\n" for i, row in enumerate(cases)),
        encoding="utf-8", newline="\n",
    )
    source = Path(__file__).with_name("cleanroom_profile_checker_v3.cpp")
    exe = output / "cleanroom_profile_checker_v3.exe"
    subprocess.run(["clang++", "-O3", "-std=c++20", "-Wall", "-Wextra", str(source), "-o", str(exe)], check=True)
    records_path = output / "cleanroom_profile_records_v3.txt"
    completed = subprocess.run([str(exe), str(case_file), str(cert), str(records_path)], capture_output=True, text=True, check=True)
    counts = {}
    for line in completed.stderr.splitlines():
        k, ordered, folded, canonical, above = map(int, line.split())
        counts[str(k)] = {
            "ordered_profiles": ordered,
            "folded_matching_profiles": folded,
            "canonical_profiles": canonical,
            "above_incumbent": above,
        }

    records = []
    for line in records_path.read_text(encoding="utf-8").splitlines():
        fields = line.split()
        records.append({"case_id": int(fields[0]), "k": int(fields[1]),
                        "value": fields[2], "profile": list(map(int, fields[3:]))})
    clean_above = {(r["k"], tuple(r["profile"]), r["value"])
                   for r in records if int(r["value"]) >= M}
    production_targets = {
        (int(r["target"]["k"]), tuple(map(int, r["target"]["correlations"])),
         r["target"]["absolute_profile_product"])
        for r in load(cert / "lift_targets.json")
    }
    missing = sorted(production_targets - clean_above)
    extra = sorted(clean_above - production_targets)
    pruned = [r for r in records if int(r["value"]) < M]
    with (output / "pruned_profiles_v3.jsonl").open("w", encoding="utf-8", newline="\n") as stream:
        for row in pruned:
            stream.write(json.dumps({
                "profile_id": digest([row["k"], row["profile"]]),
                "profile_data": {"k": row["k"], "values": row["profile"]},
                "bound_type": "exact clean-room finite-field profile product",
                "exact_upper_bound": row["value"], "incumbent": str(M),
                "pruned": True, "equivalence": "one affine correlation-profile class",
            }, sort_keys=True) + "\n")
    report = {
        "status": "PASS" if not missing and not extra else "FAIL",
        "implementation": str(source.as_posix()),
        "production_executables_called": False, "production_verifier_called": False,
        "counts_by_weight": counts, "clean_canonical_records": len(records),
        "clean_above_incumbent": len(clean_above), "production_target_count": len(production_targets),
        "missing_from_clean": [list(x) for x in missing], "extra_from_clean": [list(x) for x in extra],
        "clean_target_set_sha256": digest(sorted(clean_above)),
        "production_target_set_sha256": digest(sorted(production_targets)),
        "pruned_canonical_profile_count": len(pruned),
        "pruning_scope": "all canonical profile classes reaching the exact profile-product stage",
    }
    (output / "cleanroom_profile_report_v3.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("CLEAN-ROOM PROFILE EVALUATOR")
    for k in sorted(counts): print(k, counts[k])
    print("missing_from_clean =", len(missing)); print("extra_from_clean =", len(extra))
    print("clean_above_incumbent =", len(clean_above)); print("production_target_count =", len(production_targets))
    print("INDEPENDENT CORRELATION PROFILE CHECK:", report["status"])
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=Path("certificates/order55_global"))
    parser.add_argument("--output", type=Path, default=Path("certificates/order55_global_v2"))
    args = parser.parse_args()
    raise SystemExit(run(args.certificate.resolve(), args.output.resolve()))
