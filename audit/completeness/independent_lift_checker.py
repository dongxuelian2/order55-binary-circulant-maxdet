"""Independent structural audit of the order-55 lift partition.

This checker validates the mathematical task objects, their disjoint file
partition, and the recorded witness payloads.  It deliberately does not call
the production verifier or the native lift executable.  It therefore reports
the stronger global fiber-union statement separately as NOT_COMPLETED.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

P = 2305843009213696591
Q = 2305843009213697141
N = 55
M = 134694094094758395331307111329132


def digest(value) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def correlations(word: str) -> tuple[int, ...]:
    bits = tuple(int(bit) for bit in word)
    return tuple(sum(bits[t] * bits[(t + s) % N] for t in range(N)) for s in range(N))


def mul(a: int, b: int, p: int) -> int:
    return a * b % p


def primitive_root(p: int) -> int:
    for candidate in range(2, p):
        root = pow(candidate, (p - 1) // N, p)
        if pow(root, 5, p) != 1 and pow(root, 11, p) != 1:
            return root
    raise RuntimeError("no primitive 55th root found")


def determinant_mod(word: str, p: int, root: int) -> int:
    bits = tuple(int(bit) for bit in word)
    result = 1
    for j in range(N):
        eigenvalue = sum(bit * pow(root, j * t, p) for t, bit in enumerate(bits)) % p
        result = mul(result, eigenvalue, p)
    return result


def absolute_determinant(word: str, roots: tuple[int, int]) -> int:
    residues = (determinant_mod(word, P, roots[0]), determinant_mod(word, Q, roots[1]))
    delta = (residues[1] - residues[0]) * pow(P, -1, Q) % Q
    value = residues[0] + P * delta
    if value > P * Q // 2:
        value -= P * Q
    return abs(value)


def task_payload(target: dict, row: tuple[int, ...], col: tuple[int, ...]) -> tuple[int, ...]:
    corr = tuple(target["correlations"])
    return tuple(row) + tuple(col) + (target["k"],) + corr + tuple(reversed(corr))


def parse_task(line: str) -> tuple[str, tuple[int, ...]]:
    fields = line.split()
    if len(fields) != 72:
        raise ValueError(f"expected 72 fields, got {len(fields)}")
    return fields[0], tuple(map(int, fields[1:]))


def run(certificate: Path, output: Path) -> int:
    output.mkdir(parents=True, exist_ok=True)
    coverage = load(certificate / "lift_targets.json")
    expected: dict[str, tuple[int, ...]] = {}
    expected_target: dict[str, int] = {}
    target_task_counts = {}
    for target_id, entry in enumerate(coverage):
        rows = [tuple(row) for row in entry["row_profiles"]]
        cols = [tuple(col) for col in entry["column_profiles"]]
        target_task_counts[str(target_id)] = len(rows) * len(cols)
        for row_index, row in enumerate(rows):
            for col_index, col in enumerate(cols):
                task = f"{target_id}_{row_index}_{col_index}"
                expected[task] = task_payload(entry["target"], row, col)
                expected_target[task] = target_id

    actual: dict[str, tuple[int, ...]] = {}
    file_for_task: dict[str, str] = {}
    duplicate_tasks = []
    malformed = []
    for path in sorted(certificate.glob("lift_part[0-9].txt")):
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            try:
                task, payload = parse_task(line)
            except ValueError as exc:
                malformed.append({"file": path.name, "line": line_number, "error": str(exc)})
                continue
            if task in actual:
                duplicate_tasks.append(task)
            actual[task] = payload
            file_for_task[task] = path.name

    expected_set = set(expected)
    actual_set = set(actual)
    missing = sorted(expected_set - actual_set)
    extra = sorted(actual_set - expected_set)
    payload_mismatches = sorted(task for task in expected_set & actual_set if actual[task] != expected[task])

    audit_task_sets = {}
    audit_mismatches = []
    audit_rows = 0
    for path in sorted(certificate.glob("lift_part[0-9]_audit.jsonl")):
        seen = []
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            row = json.loads(line)
            task = row["task"]
            seen.append(task)
            audit_rows += 1
            if task not in actual or file_for_task.get(task) != path.name.replace("_audit.jsonl", ".txt"):
                audit_mismatches.append({"file": path.name, "line": line_number, "task": task})
            if row["solutions"] < 0 or row["joined_words"] < 0:
                audit_mismatches.append({"file": path.name, "line": line_number, "task": task, "reason": "negative counter"})
        audit_task_sets[path.name] = sorted(seen)
    audit_expected = set(actual)
    audit_actual = {task for tasks in audit_task_sets.values() for task in tasks}
    audit_missing = sorted(audit_expected - audit_actual)
    audit_extra = sorted(audit_actual - audit_expected)

    witnesses = []
    witness_errors = []
    roots = (primitive_root(P), primitive_root(Q))
    for path in sorted(certificate.glob("lift_part[0-9]_words.txt")):
        for line_number, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            if not line.strip():
                continue
            fields = line.split()
            if len(fields) != 2:
                witness_errors.append({"file": path.name, "line": line_number, "reason": "malformed witness"})
                continue
            task, word = fields
            row = {"file": path.name, "line": line_number, "task": task, "word": word}
            if task not in expected:
                witness_errors.append({**row, "reason": "unknown task"})
                continue
            target_id = expected_target[task]
            target = coverage[target_id]["target"]
            corr = correlations(word)
            complement = "".join("1" if bit == "0" else "0" for bit in word)
            det = absolute_determinant(complement, roots)
            row.update({"weight": word.count("1"), "correlations": list(corr[1:28]),
                       "determinant_of_complement": str(det)})
            if len(word) != N or any(bit not in "01" for bit in word):
                witness_errors.append({**row, "reason": "not a binary order-55 word"})
            elif list(corr[1:28]) != target["correlations"]:
                witness_errors.append({**row, "reason": "correlation mismatch"})
            elif det != int(target["absolute_profile_product"]):
                witness_errors.append({**row, "reason": "determinant mismatch"})
            witnesses.append(row)

    task_counts_by_part = Counter(file_for_task.values())
    target_counts_match = all(target_task_counts[str(i)] == sum(1 for task in actual if expected_target.get(task) == i)
                              for i in range(len(coverage)))
    listed_partition_pass = (
        len(coverage) == 2316 and len(expected) == 16084 and len(actual) == len(expected)
        and not missing and not extra and not duplicate_tasks and not malformed
        and not payload_mismatches and target_counts_match
        and not audit_missing and not audit_extra and not audit_mismatches
        and not witness_errors
    )
    report = {
        "status": "PASS" if listed_partition_pass else "FAIL",
        "implementation": "audit/completeness/independent_lift_checker.py",
        "production_executables_called": False,
        "production_verifier_called": False,
        "targets": len(coverage),
        "expected_tasks": len(expected),
        "actual_tasks": len(actual),
        "task_set_sha256": digest(sorted((task, list(payload)) for task, payload in actual.items())),
        "expected_task_set_sha256": digest(sorted((task, list(payload)) for task, payload in expected.items())),
        "missing_tasks": missing,
        "extra_tasks": extra,
        "duplicate_tasks": sorted(set(duplicate_tasks)),
        "malformed_tasks": malformed,
        "payload_mismatches": payload_mismatches,
        "audit_rows": audit_rows,
        "audit_missing_tasks": audit_missing,
        "audit_extra_tasks": audit_extra,
        "audit_mismatches": audit_mismatches,
        "task_counts_by_file": dict(sorted(task_counts_by_part.items())),
        "target_task_counts_match": target_counts_match,
        "witness_count": len(witnesses),
        "witness_errors": witness_errors,
        "witnesses": witnesses,
        "listed_task_partition": "PASS" if listed_partition_pass else "FAIL",
        "same_executable_replay": "split-dependent cross-check only",
        "independent_full_fiber_union": "NOT_COMPLETED",
        "global_completeness": "INCOMPLETE",
        "maximum": str(M),
    }
    report_path = output / "independent_lift_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("LIFT TASK PARTITION")
    print("targets =", len(coverage))
    print("expected_tasks =", len(expected))
    print("actual_tasks =", len(actual))
    print("audit_rows =", audit_rows)
    print("witnesses =", len(witnesses))
    print("LISTED TASK PARTITION:", report["listed_task_partition"])
    print("INDEPENDENT FULL FIBER UNION:", report["independent_full_fiber_union"])
    return 0 if listed_partition_pass else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=Path("certificates/order55_global"))
    parser.add_argument("--output", type=Path, default=Path("certificates/order55_global_v2"))
    args = parser.parse_args()
    raise SystemExit(run(args.certificate.resolve(), args.output.resolve()))
