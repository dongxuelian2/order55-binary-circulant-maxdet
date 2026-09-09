"""Small composite-order regression for the order-55 proof architecture.

The ground truth is a clean-room C++ Gray-code enumeration of every word at
orders 15 and 21 with exact two-prime CRT determinants.  This Python layer
independently checks the complement map, correlation/fold identities, formal
profile containment, affine symmetry of every ground-truth maximizer, and a
direct column-by-column reconstruction of every distinct maximizing fibre.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from collections import defaultdict
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GROUND = Path(__file__).with_name("small_order_ground_truth.cpp")
PRIMES = (1000000891, 1000001311)


def half_correlation(mask: int, order: int) -> tuple[int, ...]:
    full_mask = (1 << order) - 1
    return tuple((mask & (((mask << shift) | (mask >> (order - shift))) & full_mask)).bit_count()
                 for shift in range(1, (order + 1) // 2))


def full_correlation(mask: int, order: int) -> tuple[int, ...]:
    full_mask = (1 << order) - 1
    return (mask.bit_count(),) + tuple(
        (mask & (((mask << shift) | (mask >> (order - shift))) & full_mask)).bit_count()
        for shift in range(1, order)
    )


def fold(mask: int, order: int, modulus: int) -> tuple[int, ...]:
    return tuple(sum((mask >> index) & 1 for index in range(residue, order, modulus))
                 for residue in range(modulus))


def folded_correlation(vector: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sum(vector[index] * vector[(index + shift) % len(vector)]
                     for index in range(len(vector)))
                 for shift in range(len(vector)))


def affine_orbit(word: str, order: int) -> set[str]:
    return {
        "".join(word[(unit * index + translation) % order]
                for index in range(order))
        for unit in range(1, order)
        if gcd(unit, order) == 1
        for translation in range(order)
    }


def bareiss(matrix: list[list[int]]) -> int:
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(len(a) - 1):
        pivot_row = next((row for row in range(pivot_index, len(a))
                          if a[row][pivot_index]), None)
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            a[pivot_index], a[pivot_row] = a[pivot_row], a[pivot_index]
            sign = -sign
        pivot = a[pivot_index][pivot_index]
        for row in range(pivot_index + 1, len(a)):
            for column in range(pivot_index + 1, len(a)):
                numerator = (a[row][column] * pivot
                             - a[row][pivot_index] * a[pivot_index][column])
                if numerator % previous:
                    raise AssertionError("non-exact Bareiss division")
                a[row][column] = numerator // previous
        for row in range(pivot_index + 1, len(a)):
            a[row][pivot_index] = 0
        previous = pivot
    return sign * a[-1][-1]


def matrix(word: str) -> list[list[int]]:
    order = len(word)
    return [[int(word[(column - row) % order]) for column in range(order)]
            for row in range(order)]


def formal_profiles(order: int, weight: int) -> set[tuple[int, ...]]:
    half = (order - 1) // 2
    total = weight * (weight - 1) // 2
    result: set[tuple[int, ...]] = set()
    values = [0] * half

    def visit(position: int, minimum: int, remaining: int):
        if position == half:
            if remaining == 0:
                result.add(tuple(values))
            return
        slots = half - position
        if remaining < minimum * slots:
            return
        for value in range(minimum, min(weight, remaining // slots) + 1):
            values[position] = value
            visit(position + 1, value, remaining - value)

    visit(0, 0, total)
    return result


def crt(first: int, second: int) -> int:
    p, q = PRIMES
    value = first + p * (((second - first) * pow(p, -1, q)) % q)
    if value > p * q // 2:
        value -= p * q
    return value


def run_ground_truth(order: int) -> dict:
    with tempfile.TemporaryDirectory(prefix="small-order-audit-") as directory:
        executable = Path(directory) / "small_order_ground_truth.exe"
        subprocess.run(["clang++", "-O2", "-std=c++20", str(GROUND), "-o", str(executable)],
                       check=True)
        raw = subprocess.check_output([str(executable), str(order)], text=True)
    return json.loads(raw)


def torus_positions(order: int, rows: int, columns: int) -> dict[tuple[int, int], int]:
    row_coefficient = columns * pow(columns, -1, rows) % order
    column_coefficient = rows * pow(rows, -1, columns) % order
    positions = {
        (row, column): (row_coefficient * row + column_coefficient * column) % order
        for row in range(rows) for column in range(columns)
    }
    if len(set(positions.values())) != order:
        raise AssertionError("CRT torus coordinates are not bijective")
    return positions


def reconstruct_fibre(order: int, rows: int, columns: int,
                      row_margin: tuple[int, ...],
                      column_margin: tuple[int, ...],
                      target: tuple[int, ...]) -> int:
    positions = torus_positions(order, rows, columns)
    choices = {
        column: [mask for mask in range(1 << rows)
                 if mask.bit_count() == column_margin[column]]
        for column in range(columns)
    }
    partial = [0] * rows
    count = 0

    def visit(column: int, word: int):
        nonlocal count
        if column == columns:
            if tuple(partial) != row_margin:
                return
            if half_correlation(word, order) == target:
                count += 1
            return
        for mask in choices[column]:
            updated = list(partial)
            for row in range(rows):
                updated[row] += (mask >> row) & 1
            if any(updated[row] > row_margin[row] for row in range(rows)):
                continue
            partial[:] = updated
            next_word = word
            for row in range(rows):
                if (mask >> row) & 1:
                    next_word |= 1 << positions[row, column]
            visit(column + 1, next_word)
            partial[:] = [partial[row] - ((mask >> row) & 1)
                          for row in range(rows)]

    visit(0, 0)
    return count


def audit_order(order: int, ground: dict) -> dict:
    factors = (3, order // 3)
    total = 1 << order
    actual_profiles: dict[int, set[tuple[int, ...]]] = defaultdict(set)
    max_keys = set()
    identity_checks = 0
    for mask in range(total):
        weight = mask.bit_count()
        lower_mask = mask if weight <= order // 2 else mask ^ (total - 1)
        lower_weight = lower_mask.bit_count()
        half = half_correlation(lower_mask, order)
        full = full_correlation(lower_mask, order)
        if 2 * sum(half) != lower_weight * (lower_weight - 1):
            raise AssertionError("correlation first moment failed")
        actual_profiles[lower_weight].add(tuple(sorted(half)))
        for modulus in factors:
            margin = fold(lower_mask, order, modulus)
            aggregate = tuple(sum(full[shift] for shift in range(order)
                                  if shift % modulus == residue)
                              for residue in range(modulus))
            if aggregate != folded_correlation(margin):
                raise AssertionError("folded correlation identity failed")
        identity_checks += 1

    for weight, profiles in actual_profiles.items():
        if not profiles <= formal_profiles(order, weight):
            raise AssertionError(("actual profile missing from formal universe", order, weight))

    maximizers = set(ground["maximizers"])
    for word in maximizers:
        exact = bareiss(matrix(word))
        if exact != ground["maximum"]:
            raise AssertionError(("small ground truth matrix check failed", order, word, exact))
        weight = word.count("1")
        complement = "".join("1" if bit == "0" else "0" for bit in word)
        complement_det = bareiss(matrix(complement))
        if weight * complement_det != (order - weight) * exact:
            raise AssertionError("small complement identity failed")
        orbit = affine_orbit(word, order)
        if not orbit <= maximizers:
            raise AssertionError("affine image of a maximizer is not a maximizer")
        # The string's character zero is coordinate zero in the C++ output;
        # construct the integer without reversing to preserve that convention.
        mask = sum((bit == "1") << index for index, bit in enumerate(word))
        half = half_correlation(mask, order)
        row_margin = fold(mask, order, factors[0])
        column_margin = fold(mask, order, factors[1])
        max_keys.add((half, row_margin, column_margin))

    fibres = []
    for half, row_margin, column_margin in sorted(max_keys):
        fibres.append({
            "row_margin": list(row_margin),
            "column_margin": list(column_margin),
            "solutions": reconstruct_fibre(order, factors[0], factors[1],
                                            row_margin, column_margin, half),
        })
    if any(fibre["solutions"] == 0 for fibre in fibres):
        raise AssertionError("maximizing fibre was not reconstructed")

    weight_best = ground["weight_best"]
    for weight in range(1, order // 2 + 1):
        left = int(weight_best[weight]["maximum"])
        right = int(weight_best[order - weight]["maximum"])
        if weight * right != (order - weight) * left:
            raise AssertionError(("weight-stratum complement maxima failed", order, weight))
    return {
        "order": order,
        "enumerated_words": total,
        "maximum": ground["maximum"],
        "maximizers": len(maximizers),
        "affine_orbit_checks": len(maximizers),
        "identity_checks": identity_checks,
        "actual_profile_counts": {str(k): len(v) for k, v in actual_profiles.items()},
        "formal_profile_counts": {str(k): len(formal_profiles(order, k))
                                   for k in actual_profiles},
        "maximizing_fibres": fibres,
        "status": "PASS",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("small_order_audit.json"))
    args = parser.parse_args()
    reports = {}
    for order in (15, 21):
        ground = run_ground_truth(order)
        reports[str(order)] = audit_order(order, ground)
    args.output.write_text(json.dumps(reports, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "orders": [15, 21],
        "status": "PASS",
        "maxima": {order: reports[str(order)]["maximum"] for order in (15, 21)},
        "words": {order: reports[str(order)]["enumerated_words"] for order in (15, 21)},
    }, sort_keys=True))


if __name__ == "__main__":
    main()
