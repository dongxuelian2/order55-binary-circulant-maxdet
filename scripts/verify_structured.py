"""Independently verify the order-55 fixed-weight circulant certificate."""

from __future__ import annotations

import argparse
from math import comb
import json
from pathlib import Path

import sympy as sp


def circulant(word: str) -> sp.Matrix:
    values = [int(value) for value in word]
    order = len(values)
    return sp.Matrix(order, order, lambda row, column: values[(column - row) % order])


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "certificate",
        nargs="?",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "certificates" / "order55_fixed_weight_1_7.json",
    )
    args = parser.parse_args()
    data = json.loads(args.certificate.read_text(encoding="utf-8"))
    order = int(data["order"])
    enumeration = data["enumeration"]
    prime1 = int(enumeration["prime_1"])
    prime2 = int(enumeration["prime_2"])
    assert order == 55
    assert sp.isprime(prime1) and sp.isprime(prime2)
    assert (prime1 - 1) % order == 0 and (prime2 - 1) % order == 0

    total = 0
    for result in data["results"]:
        weight = int(result["weight"])
        word = result["word"]
        expected = int(result["absolute_determinant"])
        representatives = int(result["representatives_evaluated"])
        assert len(word) == order and set(word) <= {"0", "1"}
        assert word.count("1") == weight and word[0] == "1"
        assert representatives == comb(order - 1, weight - 1)
        total += representatives
        determinant = int(circulant(word).det(method="domain-ge"))
        assert abs(determinant) == expected
        assert determinant % prime1 == int(result["residue_prime_1"])
        assert determinant % prime2 == int(result["residue_prime_2"])
        assert expected * expected <= weight**order
        assert (prime1 * prime2 // 2) ** 2 > weight**order
        print(f"weight={weight} determinant={expected} representatives={representatives} PASS")

    assert total == int(enumeration["total_representatives"])
    print(f"order={order} total_representatives={total} PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
