"""Independent affine-orbit audit; complements need not be canonical representatives."""

from __future__ import annotations

import argparse
import json
from math import gcd
from pathlib import Path

N = 55
P = 2305843009213696591
Q = 2305843009213697141
EXPECTED_WORD = "0000000100011111011011001101011011010100011110101000111"
EXPECTED_COMPLEMENT = "1111111011100000100100110010100100101011100001010111000"
EXPECTED_MAX = 134694094094758395331307111329132
EXPECTED_COMPLEMENT_DET = 129883590734231309783760428781663


def load(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def affine_images(word: str) -> set[str]:
    return {"".join(word[(u * i + t) % N] for i in range(N))
            for u in range(N) if gcd(u, N) == 1 for t in range(N)}


def canonical(word: str) -> str:
    return min(affine_images(word))


def mul(a: int, b: int, p: int) -> int:
    return a * b % p


def root55(p: int) -> int:
    for candidate in range(2, p):
        root = pow(candidate, (p - 1) // N, p)
        if pow(root, 5, p) != 1 and pow(root, 11, p) != 1:
            return root
    raise RuntimeError("primitive 55th root not found")


def determinant_mod(word: str, p: int, root: int) -> int:
    result = 1
    for j in range(N):
        value = sum(int(bit) * pow(root, j * t, p) for t, bit in enumerate(word)) % p
        result = mul(result, value, p)
    return result


def absolute_determinant(word: str, roots: tuple[int, int]) -> int:
    r = determinant_mod(word, P, roots[0])
    s = determinant_mod(word, Q, roots[1])
    value = r + P * ((s - r) * pow(P, -1, Q) % Q)
    if value > P * Q // 2:
        value -= P * Q
    return abs(value)


def run(certificate: Path, output: Path) -> int:
    recorded = load(certificate / "winner.json")["classes"][0]
    word, complement = recorded["word"], recorded["complement"]
    orbit, complement_orbit = affine_images(word), affine_images(complement)
    lifted = load(certificate / "lifted_words.json")
    lifted_canonicals = sorted({canonical(row["word"]) for row in lifted})
    roots = (root55(P), root55(Q))
    winner_det = absolute_determinant(word, roots)
    complement_det = absolute_determinant(complement, roots)
    passed = (
        word == EXPECTED_WORD and complement == EXPECTED_COMPLEMENT
        and len(orbit) == 2200 and recorded["orbit_size"] == 2200
        and recorded["stabilizer_size"] == 1 and canonical(word) == word
        and len(complement_orbit) == 2200 and lifted_canonicals == [word]
        and winner_det == EXPECTED_MAX and complement_det == EXPECTED_COMPLEMENT_DET
    )
    report = {
        "status": "PASS" if passed else "FAIL",
        "implementation": "audit/completeness/affine_canonicalization_audit_v2.py",
        "production_executables_called": False, "production_verifier_called": False,
        "winner_word": word, "complement": complement,
        "winner_canonical": canonical(word), "complement_canonical": canonical(complement),
        "winner_orbit_size_independent": len(orbit),
        "complement_orbit_size_independent": len(complement_orbit),
        "winner_stabilizer_size_independent": 2200 // len(orbit),
        "recorded_orbit_size": recorded["orbit_size"],
        "recorded_stabilizer_size": recorded["stabilizer_size"],
        "lifted_witness_canonical_classes": lifted_canonicals,
        "winner_determinant_independent": str(winner_det),
        "complement_determinant_independent": str(complement_det),
        "winner_weight": word.count("1"), "complement_weight": complement.count("1"),
        "complement_is_not_required_to_be_canonical": True,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "affine_canonicalization_report_v2.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print("AFFINE CANONICALIZATION")
    print("winner_orbit_size =", len(orbit)); print("winner_stabilizer_size =", 2200 // len(orbit))
    print("complement_orbit_size =", len(complement_orbit)); print("lifted_witness_classes =", len(lifted_canonicals))
    print("winner_determinant =", winner_det); print("complement_determinant =", complement_det)
    print("INDEPENDENT AFFINE AUDIT:", report["status"])
    return 0 if passed else 1


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--certificate", type=Path, default=Path("certificates/order55_global"))
    parser.add_argument("--output", type=Path, default=Path("certificates/order55_global_v2"))
    args = parser.parse_args()
    raise SystemExit(run(args.certificate.resolve(), args.output.resolve()))
