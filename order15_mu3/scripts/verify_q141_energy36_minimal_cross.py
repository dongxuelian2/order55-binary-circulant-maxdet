"""Exact Schur enumeration for the two tight Q=141 energy-36 supports."""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from itertools import product

from maxdet.mu3 import Eisenstein, determinant, inner_product_values

try:
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_trace_stability import stationary_upper


def _minor(matrix, removed_row: int, removed_column: int):
    return [
        [value for column, value in enumerate(row) if column != removed_column]
        for row_index, row in enumerate(matrix)
        if row_index != removed_row
    ]


def _adjugate(matrix: list[list[Eisenstein]]) -> list[list[Eisenstein]]:
    size = len(matrix)
    return [
        [
            determinant(_minor(matrix, column, row)) * (-1 if (row + column) % 2 else 1)
            for column in range(size)
        ]
        for row in range(size)
    ]


def _quadratic(left, matrix, right) -> Eisenstein:
    return sum(
        (left[i].conjugate() * matrix[i][j] * right[j] for i in range(len(left)) for j in range(len(right))),
        Eisenstein(),
    )


def _isolate_sums(values: tuple[Eisenstein, ...], isolates: int) -> tuple[Eisenstein, ...]:
    products = sorted(
        {left.conjugate() * right for left in values for right in values},
        key=lambda value: (value.a, value.b),
    )
    assert len(products) == 3
    return tuple(
        first * products[0] + second * products[1] + (isolates - first - second) * products[2]
        for first in range(isolates + 1)
        for second in range(isolates + 1 - first)
    )


def _profile_maximum(active_vertices: int, tree_edges, cycle_edge, isolates: int) -> dict[str, int]:
    catalogue = inner_product_values(15)

    def delta(entry) -> int:
        _a, b, c = entry["counts"][0]
        return (b + 2 * c) % 3

    cross_values = tuple(
        entry["value"] for entry in catalogue if entry["norm"] == 3 and delta(entry) == 1
    )
    cycle_values = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    outside_values = tuple(entry["value"] for entry in catalogue if entry["norm"] == 27)
    vectors = tuple(product(cross_values, repeat=active_vertices))
    isolate_sums = _isolate_sums(cross_values, isolates)
    maximum = Fraction(0)
    distinct_cross_terms = set()
    for cycle_value in cycle_values:
        block = [
            [Eisenstein(15 if i == j else 0) for j in range(active_vertices)]
            for i in range(active_vertices)
        ]
        for i, j in tree_edges:
            block[i][j] = block[j][i] = Eisenstein(3)
        i, j = cycle_edge
        block[i][j] = cycle_value
        block[j][i] = cycle_value.conjugate()
        block_determinant = determinant(block)
        assert block_determinant.b == 0 and block_determinant.a > 0
        determinant_value = block_determinant.a
        adjugate = _adjugate(block)
        diagonals = []
        for vector in vectors:
            numerator = _quadratic(vector, adjugate, vector)
            assert numerator.b == 0
            diagonals.append(Fraction(15) - Fraction(isolates, 5) - Fraction(numerator.a, determinant_value))
        residual_cache = {}
        for left_index, left in enumerate(vectors):
            for right_index, right in enumerate(vectors):
                cross = _quadratic(left, adjugate, right)
                key = (cross.a, cross.b)
                distinct_cross_terms.add((determinant_value, *key))
                if key not in residual_cache:
                    residual_cache[key] = min(
                        (
                            Fraction(outside.a) - Fraction(cross.a, determinant_value) - Fraction(isolate_sum.a, 15)
                        ) ** 2
                        - (
                            Fraction(outside.a) - Fraction(cross.a, determinant_value) - Fraction(isolate_sum.a, 15)
                        ) * (
                            Fraction(outside.b) - Fraction(cross.b, determinant_value) - Fraction(isolate_sum.b, 15)
                        )
                        + (
                            Fraction(outside.b) - Fraction(cross.b, determinant_value) - Fraction(isolate_sum.b, 15)
                        ) ** 2
                        for outside in outside_values
                        for isolate_sum in isolate_sums
                    )
                schur = diagonals[left_index] * diagonals[right_index] - residual_cache[key]
                candidate = Fraction(determinant_value * 15**isolates) * schur
                if candidate > maximum:
                    maximum = candidate
    assert maximum.denominator == 1
    return {
        "maximum": maximum.numerator,
        "phase_assignments": len(cycle_values) * len(vectors) ** 2 * len(outside_values) * len(isolate_sums),
        "distinct_cross_terms": len(distinct_cross_terms),
    }


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    profiles = {
        "K2+K3": _profile_maximum(
            5, ((0, 1), (2, 3), (2, 4)), (3, 4), 8
        ),
        "paw": _profile_maximum(
            4, ((0, 1), (0, 2), (0, 3)), (1, 2), 9
        ),
    }
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert profiles["K2+K3"]["maximum"] == 248427409860000000
    assert all(profile["maximum"] < q153 for profile in profiles.values())
    return {
        "profiles": profiles,
        "maximum": max(profile["maximum"] for profile in profiles.values()),
        "over_q153": float(max(profile["maximum"] for profile in profiles.values()) / q153),
        "theorem": "Both tight Q=141 energy-36 minimal-cross supports lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
