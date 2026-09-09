"""Run the clean-room typed-key lift solver over all committed lift jobs."""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import tempfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
CERT = ROOT / "certificates" / "order55_global"
SOURCE = Path(__file__).with_name("independent_torus_lift.cpp")


def read_jsonl(path: Path):
    return [json.loads(line) for line in path.read_text().splitlines() if line]


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_part(executable: Path, part: int, temp: Path):
    input_path = CERT / f"lift_part{part}.txt"
    audit_path = temp / f"lift_part{part}_audit.jsonl"
    words_path = temp / f"lift_part{part}_words.txt"
    subprocess.run([str(executable), str(input_path), str(audit_path),
                    str(words_path), "6"], cwd=ROOT, check=True)
    new_audits = read_jsonl(audit_path)
    old_audits = read_jsonl(CERT / f"lift_part{part}_audit.jsonl")
    old_by_task = {row["task"]: row for row in old_audits}
    new_by_task = {row["task"]: row for row in new_audits}
    if set(new_by_task) != set(old_by_task):
        raise AssertionError(("task-id set mismatch", part))
    for task, row in new_by_task.items():
        old = old_by_task[task]
        for key in ("joined_words", "solutions"):
            if row[key] != old[key]:
                raise AssertionError(("lift count mismatch", part, task, key,
                                      row[key], old[key]))
    new_words = set(words_path.read_text().splitlines())
    old_words = set((CERT / f"lift_part{part}_words.txt").read_text().splitlines())
    if new_words != old_words:
        raise AssertionError(("solution set mismatch", part,
                              len(new_words), len(old_words)))
    return {
        "part": part,
        "tasks": len(new_by_task),
        "joined_words": sum(row["joined_words"] for row in new_audits),
        "solutions": sum(row["solutions"] for row in new_audits),
        "audit_sha256": sha256(audit_path),
        "words_sha256": sha256(words_path),
        "split": 6,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("independent_lift_audit.json"))
    args = parser.parse_args()
    if not 1 <= args.workers <= 8:
        raise SystemExit("workers must be in 1..8")
    with tempfile.TemporaryDirectory(prefix="order55-clean-lift-") as directory:
        temp = Path(directory)
        executable = temp / "independent_torus_lift.exe"
        subprocess.run(["clang++", "-O3", "-std=c++20", "-Wall", "-Wextra",
                        "-Wconversion", str(SOURCE), "-o", str(executable)],
                       cwd=ROOT, check=True)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            summaries = list(pool.map(lambda part: run_part(executable, part, temp), range(8)))
    summaries.sort(key=lambda row: row["part"])
    output = {
        "status": "PASS",
        "solver": "audit/order55/independent_torus_lift.cpp",
        "split": 6,
        "parts": summaries,
        "tasks": sum(row["tasks"] for row in summaries),
        "joined_words": sum(row["joined_words"] for row in summaries),
        "solutions": sum(row["solutions"] for row in summaries),
    }
    args.output.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "status": "PASS",
        "tasks": output["tasks"],
        "joined_words": output["joined_words"],
        "solutions": output["solutions"],
        "split": output["split"],
    }, sort_keys=True))


if __name__ == "__main__":
    main()
