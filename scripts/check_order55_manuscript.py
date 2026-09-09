"""Check publication text against the authoritative order-55 certificate.

This is a packaging-time transcription check.  It does not prove the global
theorem; it prevents the publication source from silently diverging from the
already-frozen certificate.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path


def required(text: str, needle: str, label: str) -> None:
    if needle not in text:
        raise SystemExit(f"MISSING {label}: {needle}")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("manuscript", type=Path)
    parser.add_argument(
        "--certificate",
        type=Path,
        default=Path("certificates/order55_global"),
    )
    args = parser.parse_args()

    manuscript = args.manuscript.read_text(encoding="utf-8")
    folder = args.certificate
    global_record = json.loads((folder / "global.json").read_text(encoding="utf-8"))
    winner_record = json.loads((folder / "winner.json").read_text(encoding="utf-8"))
    winner = winner_record["classes"][0]

    if global_record["status"] != "COMPLETE":
        raise SystemExit("certificate is not marked COMPLETE")
    if len(winner["word"]) != 55 or set(winner["word"]) - {"0", "1"}:
        raise SystemExit("certificate winner is not a 55-bit binary word")

    complement = "".join("1" if bit == "0" else "0" for bit in winner["word"])
    expected = {
        "maximum": str(global_record["maximum"]),
        "winner": winner["word"],
        "weight": str(winner["weight"]),
        "complement": complement,
        "complement_determinant": "129883590734231309783760428781663",
        "maximizing_words": str(global_record["maximizing_words"]),
        "affine_class_count": str(global_record["affine_class_count"]),
        "stabilizer_size": str(winner["stabilizer_size"]),
        "correlation_profiles": str(global_record["correlation_targets"]),
        "explicit_profiles": "6845230",
        "lift_tasks": str(global_record["lift_tasks"]),
        "joined_words": str(global_record["joined_words"]),
    }
    labels = {
        "maximum": "maximum",
        "winner": "canonical winner",
        "weight": "winner weight",
        "complement": "complement",
        "complement_determinant": "complement determinant",
        "maximizing_words": "maximizing words",
        "affine_class_count": "affine class count",
        "stabilizer_size": "stabilizer size",
        "correlation_profiles": "correlation target count",
        "explicit_profiles": "explicit profile count",
        "lift_tasks": "lift task count",
        "joined_words": "joined-word count",
    }
    for key, value in expected.items():
        required(manuscript, value, labels[key])

    required(manuscript, "Qichao Wang", "author Qichao Wang")
    required(manuscript, "Hebei University of Technology", "Qichao affiliation")
    required(manuscript, "Daoyu Dong", "author Daoyu Dong")
    required(
        manuscript,
        "University of Electronic Science and Technology of China",
        "Daoyu affiliation",
    )
    required(manuscript, "contributed equally", "equal contribution statement")
    required(manuscript, "D_{01}(55)", "D_01 notation")

    for bad in ("TODO", "TBD", "YOUR EMAIL", "INSERT HERE"):
        if re.search(re.escape(bad), manuscript, re.IGNORECASE):
            raise SystemExit(f"UNRESOLVED PLACEHOLDER: {bad}")

    print("ORDER 55 MANUSCRIPT CONSISTENCY: PASS")
    print(f"maximum = {expected['maximum']}")
    print(f"canonical winner length = {len(expected['winner'])}")
    print(f"maximizing words = {expected['maximizing_words']}")
    print(f"explicit profiles = {expected['explicit_profiles']}")
    print(f"lift tasks = {expected['lift_tasks']}")
    print(f"joined words = {expected['joined_words']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
