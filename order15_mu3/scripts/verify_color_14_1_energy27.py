"""Exact abstract-Gram enumeration for (14,1,0) at leaf energy 27."""

from __future__ import annotations

import json
from fractions import Fraction
from itertools import product

from maxdet.mu3 import Eisenstein, OMEGA, OMEGA2, ONE, determinant, inner_product_values

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm, valuation
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_bounds import is_rational_eisenstein_norm, valuation
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_trace_stability import BENCHMARK


ENERGY_LIMIT = 168
MAX_CROSS_EXCESS = 6


def arithmetic_admissible(value: int) -> bool:
    return valuation(value, 3) >= 14 and is_rational_eisenstein_norm(value)


def orbit_representatives(values: tuple[Eisenstein, ...]) -> tuple[Eisenstein, ...]:
    remaining = set(values)
    result = []
    while remaining:
        value = min(remaining, key=lambda z: (z.norm(), z.a, z.b))
        orbit = {value * root for root in (ONE, OMEGA, OMEGA2)}
        assert orbit <= set(values)
        result.append(value)
        remaining -= orbit
    return tuple(result)


def cross_vectors(values: tuple[Eisenstein, ...], length: int, max_excess: int):
    representatives = orbit_representatives(values)

    def extend(prefix, excess):
        if len(prefix) == length:
            yield tuple(prefix)
            return
        choices = representatives if not prefix else values
        for value in choices:
            increment = (value.norm() - 3) // 9
            if excess + increment <= max_excess:
                yield from extend(prefix + [value], excess + increment)

    yield from extend([], 0)


def component_types(
    vertices: int,
    edges: tuple[tuple[int, int], ...],
    edge_options: tuple[tuple[Eisenstein, ...], ...],
    cross_values: tuple[Eisenstein, ...],
) -> set[tuple[int, int, int]]:
    """Return (cross energy, det A, Schur numerator) local types."""

    result = set()
    vectors = tuple(cross_vectors(cross_values, vertices, MAX_CROSS_EXCESS))
    for labels in product(*edge_options):
        block = [[Eisenstein(15 if i == j else 0) for j in range(vertices)] for i in range(vertices)]
        for (i, j), label in zip(edges, labels):
            block[i][j] = label
            block[j][i] = label.conjugate()
        block_det = determinant(block)
        assert block_det.b == 0
        if block_det.a <= 0:
            continue
        # Compute adj(A) once for this edge labelling.  Evaluating the exact
        # quadratic form u adj(A) u* is much faster than a fresh (k+1)-order
        # Bareiss determinant for every cross vector.
        adjugate = []
        for i in range(vertices):
            adjugate_row = []
            for j in range(vertices):
                minor = [
                    [block[row][column] for column in range(vertices) if column != i]
                    for row in range(vertices) if row != j
                ]
                cofactor = determinant(minor)
                adjugate_row.append(cofactor if (i + j) % 2 == 0 else -cofactor)
            adjugate.append(adjugate_row)
        for vector in vectors:
            quadratic = sum(
                (vector[i] * adjugate[i][j] * vector[j].conjugate() for i in range(vertices) for j in range(vertices)),
                Eisenstein(),
            )
            assert quadratic.b == 0 and quadratic.a >= 0
            result.add((sum(value.norm() for value in vector), block_det.a, quadratic.a))
    return result


def reachable_sums(values: tuple[int, ...], length: int, maximum: int) -> set[int]:
    reachable = {0}
    for _ in range(length):
        reachable = {old + value for old in reachable for value in values if old + value <= maximum}
    return reachable


def evaluate_profiles(block_profiles, isolates: int, cross_norms: tuple[int, ...]) -> set[int]:
    remaining = reachable_sums(cross_norms, isolates, ENERGY_LIMIT)
    values = set()
    for h in range(MAX_CROSS_EXCESS + 1):
        cross_energy = 42 + 9 * h
        for blocks in block_profiles:
            local_energy = sum(block[0] for block in blocks)
            singles_energy = cross_energy - local_energy
            if singles_energy not in remaining:
                continue
            determinant_factor = 15**isolates
            schur = Fraction(15) - Fraction(singles_energy, 15)
            for _energy, block_det, numerator in blocks:
                determinant_factor *= block_det
                schur -= Fraction(numerator, block_det)
            value = Fraction(determinant_factor) * schur
            if value > BENCHMARK:
                assert value.denominator == 1
                values.add(value.numerator)
    return values


def main() -> None:
    catalogue = inner_product_values(15)
    cross_values = tuple(entry["value"] for entry in catalogue if 0 < entry["norm"] <= 81 and entry["norm"] % 9 == 3)
    cross_norms = tuple(sorted({value.norm() for value in cross_values}))
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    norm27 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 27)
    # Convenient representatives of the two cube-root orbits.
    sign_reps = (Eisenstein(-3), Eisenstein(3))
    norm27_reps = orbit_representatives(norm27)
    assert len(orbit_representatives(norm9)) == 2

    # The exact stationary determinant bound and lambda_max<23 show that
    # h>=7 cannot beat the record when the leaf internal energy is 27.
    assert Fraction(2 * 27 * 13, 14) < (23 - 15) ** 2
    assert stationary_upper_dimension(14, 27) * (15 - Fraction(42 + 9 * 7, 23)) < BENCHMARK

    pair9 = component_types(2, ((0, 1),), (sign_reps,), cross_values)
    pair27 = component_types(2, ((0, 1),), (norm27_reps,), cross_values)
    path3 = component_types(3, ((0, 1), (1, 2)), (sign_reps, sign_reps), cross_values)
    path4 = component_types(4, ((0, 1), (1, 2), (2, 3)), (sign_reps,) * 3, cross_values)
    star4 = component_types(4, ((0, 1), (0, 2), (0, 3)), (sign_reps,) * 3, cross_values)
    triangle3 = component_types(3, ((0, 1), (1, 2), (0, 2)), (sign_reps, sign_reps, norm9), cross_values)

    cases = {
        "one_norm27_edge": evaluate_profiles(((block,) for block in pair27), 12, cross_norms),
        "three_disjoint_norm9_edges": evaluate_profiles(product(pair9, repeat=3), 8, cross_norms),
        "path3_plus_edge": evaluate_profiles(product(path3, pair9), 9, cross_norms),
        "path4": evaluate_profiles(((block,) for block in path4), 10, cross_norms),
        "star4": evaluate_profiles(((block,) for block in star4), 10, cross_norms),
        "triangle3": evaluate_profiles(((block,) for block in triangle3), 11, cross_norms),
    }
    admissible = {name: sorted(value for value in values if arithmetic_admissible(value)) for name, values in cases.items()}
    print(json.dumps({
        "color_partition": [14, 1, 0],
        "leaf_internal_energy": 27,
        "max_cross_excess_certified": MAX_CROSS_EXCESS,
        "local_type_counts": {
            "pair_norm9": len(pair9), "pair_norm27": len(pair27), "path3": len(path3),
            "path4": len(path4), "star4": len(star4), "triangle3": len(triangle3),
        },
        "above_record_determinant_counts": {name: len(values) for name, values in cases.items()},
        "arithmetic_survivor_counts": {name: len(values) for name, values in admissible.items()},
        "distinct_arithmetic_survivors": sorted(set().union(*map(set, admissible.values()))),
    }, indent=2))


if __name__ == "__main__":
    main()
