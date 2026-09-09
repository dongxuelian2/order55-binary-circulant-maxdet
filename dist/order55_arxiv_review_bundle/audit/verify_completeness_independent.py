"""Aggregate the independently checkable order-55 completeness evidence.

This checker is deliberately read-only with respect to the production search:
it imports no production module, launches no production executable, and does
not treat a second split of the native lift as an independent generator.  It
validates the clean-room reports, the machine-readable pruning table, and the
release-gate statuses.  The current expected result is INCOMPLETE because the
full order-55 fiber-union certificate is still open.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
from pathlib import Path


INCUMBENT = 134694094094758395331307111329132


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def digest(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def require(condition: bool, message: str, failures: list[str]) -> None:
    if not condition:
        failures.append(message)


def check_pruning_table(path: Path, expected_count: int, failures: list[str]) -> dict:
    count = 0
    duplicate_ids = 0
    bad_rows = 0
    seen: set[str] = set()
    with gzip.open(path, "rt", encoding="utf-8", newline="") as stream:
        for line in stream:
            if not line.strip():
                continue
            count += 1
            try:
                row = json.loads(line)
                profile_id = row["profile_id"]
                profile_data = row["profile_data"]
                values = profile_data["values"]
                k = profile_data["k"]
                exact_upper_bound = row["exact_upper_bound"]
                incumbent = row["incumbent"]
                valid = (
                    isinstance(profile_id, str)
                    and digest([k, values]) == profile_id
                    and isinstance(values, list)
                    and len(values) == 27
                    and row["bound_type"] == "exact clean-room finite-field profile product"
                    and row["pruned"] is True
                    and int(exact_upper_bound) < INCUMBENT
                    and incumbent == str(INCUMBENT)
                )
                if profile_id in seen:
                    duplicate_ids += 1
                seen.add(profile_id)
                if not valid:
                    bad_rows += 1
            except (KeyError, TypeError, ValueError, json.JSONDecodeError):
                bad_rows += 1

    require(count == expected_count, f"pruning row count {count} != {expected_count}", failures)
    require(duplicate_ids == 0, f"duplicate pruning profile ids: {duplicate_ids}", failures)
    require(bad_rows == 0, f"invalid pruning rows: {bad_rows}", failures)
    return {
        "path": path.name,
        "rows": count,
        "expected_rows": expected_count,
        "duplicate_profile_ids": duplicate_ids,
        "invalid_rows": bad_rows,
        "schema": [
            "profile_id",
            "profile_data",
            "bound_type",
            "exact_upper_bound",
            "incumbent",
            "pruned",
        ],
        "pass": count == expected_count and duplicate_ids == 0 and bad_rows == 0,
    }


def run(audit_dir: Path, output: Path) -> int:
    failures: list[str] = []
    profile = read_json(audit_dir / "cleanroom_profile_report_v3.json")
    lift = read_json(audit_dir / "independent_lift_report.json")
    affine = read_json(audit_dir / "affine_canonicalization_report_v2.json")
    small = read_json(audit_dir / "small_order_validation_report.json")
    random = read_json(audit_dir / "random_forward_coverage_report_v6.json")
    review = read_json(audit_dir / "global_review.json")

    require(profile["status"] == "PASS", "clean-room profile comparison is not PASS", failures)
    require(profile["production_executables_called"] is False, "profile audit called a production executable", failures)
    require(profile["production_verifier_called"] is False, "profile audit called the production verifier", failures)
    require(profile["missing_from_clean"] == [], "clean-room profile set has missing records", failures)
    require(profile["extra_from_clean"] == [], "clean-room profile set has extra records", failures)
    require(profile["clean_above_incumbent"] == 2316, "clean-room target count is not 2316", failures)

    pruning = check_pruning_table(
        audit_dir / "pruned_profiles_v3.jsonl.gz",
        int(profile["pruned_canonical_profile_count"]),
        failures,
    )

    require(lift["listed_task_partition"] == "PASS", "listed lift task partition is not PASS", failures)
    require(lift["expected_tasks"] == 16084, "expected lift task count is not 16084", failures)
    require(lift["actual_tasks"] == 16084, "actual lift task count is not 16084", failures)
    require(lift["audit_rows"] == 16084, "lift audit row count is not 16084", failures)
    require(lift["production_executables_called"] is False, "lift audit called a production executable", failures)

    require(affine["status"] == "PASS", "affine audit is not PASS", failures)
    require(small["status"] == "PASS", "small-order validation is not PASS", failures)
    require(random["status"] == "PASS", "random coverage audit is not PASS", failures)
    require(random["sample_count"] == 1000000, "random coverage sample count is not 1,000,000", failures)
    require(random["unknown"] == 0, "random coverage contains UNKNOWN classifications", failures)
    require(review["release_decision"] == "MAJOR REVISION REQUIRED", "release gate is not Outcome B", failures)
    require(
        review["lift_audit"]["independent_full_fiber_union"] == "NOT_COMPLETED",
        "review gate no longer records the open independent fiber union",
        failures,
    )

    complete = not failures and review["release_decision"] == "ARXIV PACKAGE READY"
    report = {
        "status": "PASS" if complete else "INCOMPLETE" if not failures else "FAIL",
        "domain_size": review["domain_size"],
        "domain_size_meaning": "mathematical domain cardinality, not a derived coverage count",
        "production_executables_called": False,
        "checks": {
            "profile_universe": {
                "status": profile["status"],
                "clean_target_count": profile["clean_above_incumbent"],
                "missing": len(profile["missing_from_clean"]),
                "extra": len(profile["extra_from_clean"]),
            },
            "pruned_profile_table": pruning,
            "listed_lift_task_partition": lift["listed_task_partition"],
            "independent_lift_fiber_union": lift["independent_full_fiber_union"],
            "affine_canonicalization": affine["status"],
            "small_order_validation": small["status"],
            "random_forward_coverage": {
                "status": random["status"],
                "sample_count": random["sample_count"],
                "unknown": random["unknown"],
            },
        },
        "failures": failures,
        "unresolved": review["unresolved"],
    }
    output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    print("INDEPENDENT COMPLETENESS EVIDENCE")
    print("  profile_universe =", report["checks"]["profile_universe"]["status"])
    print("  pruned_profile_table =", "PASS" if pruning["pass"] else "FAIL")
    print("  listed_lift_task_partition =", lift["listed_task_partition"])
    print("  independent_full_fiber_union =", lift["independent_full_fiber_union"])
    print("  small_order_validation =", small["status"])
    print("  random_forward_coverage =", random["status"], "unknown=", random["unknown"])
    if failures:
        print("GLOBAL COMPLETENESS: FAIL")
        for failure in failures:
            print("  ", failure)
        return 1
    if complete:
        print("GLOBAL COMPLETENESS: PASS")
    else:
        print("GLOBAL COMPLETENESS: INCOMPLETE")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-dir", type=Path, default=Path("certificates/order55_global_v2"))
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("certificates/order55_global_v2/independent_completeness_report.json"),
    )
    args = parser.parse_args()
    raise SystemExit(run(args.audit_dir.resolve(), args.output.resolve()))
