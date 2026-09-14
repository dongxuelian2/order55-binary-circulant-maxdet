"""Exact rational certificate reducing counterexample row-color partitions."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, sqrt_lower
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_trace_stability import BENCHMARK, sqrt_lower


def stationary_upper_dimension(order: int, energy: int) -> Fraction:
    if energy == 0:
        return Fraction(15**order)
    candidates = []
    for multiplicity in range(1, order):
        t = sqrt_lower(Fraction(2 * energy, order * multiplicity * (order - multiplicity)))
        high = 15 + (order - multiplicity) * t
        low = 15 - multiplicity * t
        if low > 0:
            candidates.append(high**multiplicity * low ** (order - multiplicity))
    return max(candidates)


def schur_upper(size: int, internal_energy: int, lambda_upper: Fraction, cross_energy: int) -> Fraction:
    remainder = 15 - size
    trace_schur = 15 * remainder - Fraction(cross_energy, 1) / lambda_upper
    return stationary_upper_dimension(size, internal_energy) * (trace_schur / remainder) ** remainder


def partitions_of_15():
    for first in range(15, -1, -1):
        for second in range(min(first, 15 - first), -1, -1):
            third = 15 - first - second
            if 0 <= third <= second:
                yield (first, second, third)


def main() -> None:
    # If a color class of size a is internally orthogonal, Schur complement,
    # the minimum cross energy 3a(15-a), and AM-GM give this bound.
    orthogonal_class_bounds = {
        size: Fraction(15**size) * Fraction(75 - size, 5) ** (15 - size)
        for size in range(3, 13)
    }
    assert all(value < BENCHMARK for value in orthogonal_class_bounds.values())

    initial = []
    for partition in partitions_of_15():
        cross_minimum = 3 * sum(partition[i] * partition[j] for i in range(3) for j in range(i + 1, 3))
        required_internal_classes = sum(3 <= size <= 12 for size in partition)
        minimum_surviving_energy = cross_minimum + 9 * required_internal_classes
        if minimum_surviving_energy <= 168:
            initial.append((partition, cross_minimum, minimum_surviving_energy))
    assert [item[0] for item in initial] == [
        (15, 0, 0), (14, 1, 0), (13, 2, 0), (13, 1, 1),
        (12, 3, 0), (12, 2, 1), (11, 4, 0), (11, 3, 1),
        (11, 2, 2), (10, 5, 0),
    ]

    # Size 10: the budget forces one norm-9 internal edge in each class.
    assert schur_upper(10, 9, Fraction(18), 150) < BENCHMARK

    # Size 11: generic lambda_max bounds from its internal trace energy suffice.
    lambda_11 = {9: Fraction(20), 18: Fraction(21), 27: Fraction(23)}
    for energy, bound in lambda_11.items():
        assert Fraction(2 * energy * 10, 11) < (bound - 15) ** 2
        assert schur_upper(11, energy, bound, 132) < BENCHMARK

    # Size 12. At energies 9,18,27,36 the discrete same-color norm spectrum
    # gives respectively: one 3-edge; two 3-edges; one sqrt(27)-edge or three
    # 3-edges; and one 6-edge, sqrt(27)+3 edges, or four 3-edges. Entrywise
    # absolute-value comparison and the extremal small support graphs give the
    # listed lambda_max bounds. Energies 45 and 54 use the generic trace bound.
    lambda_12 = {
        9: Fraction(18),
        18: Fraction(77, 4),
        27: Fraction(21),
        36: Fraction(22),
        45: Fraction(241, 10),
        54: Fraction(25),
    }
    assert 288 < 289  # 15+3sqrt(2) < 77/4.
    assert 45 < 49    # 15+3sqrt(5) < 22 for four norm-9 edges.
    assert Fraction(165, 2) < Fraction(91, 10) ** 2
    assert 99 < 100
    for energy, bound in lambda_12.items():
        assert schur_upper(12, energy, bound, 108) < BENCHMARK

    survivors = [(15, 0, 0), (14, 1, 0), (13, 2, 0), (13, 1, 1)]
    print(json.dumps({
        "orthogonal_color_class_sizes_excluded": list(range(3, 13)),
        "energy_feasible_partitions_before_refined_schur": [item[0] for item in initial],
        "refined_largest_class_sizes_excluded": [10, 11, 12],
        "surviving_counterexample_color_partitions": survivors,
    }, indent=2))


if __name__ == "__main__":
    main()
