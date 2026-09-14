"""Exact low-internal-energy certificate for the (14,1,0) color case.

Let the 14 rows of the large color class form the leaf block A and let the
remaining row be the hub.  Cross-color Gram norms are 3 modulo 9, whereas
same-color nonzero norms are divisible by 9.  This script treats leaf
internal energy 0 and 9 exactly.
"""

from __future__ import annotations

import json
from fractions import Fraction

from maxdet.mu3 import Eisenstein, inner_product_values

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm, valuation
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_bounds import is_rational_eisenstein_norm, valuation
    from verify_trace_stability import BENCHMARK


ENERGY_LIMIT = 168


def reachable_sums(values: tuple[int, ...], length: int, maximum: int) -> set[int]:
    reachable = {0}
    for _ in range(length):
        reachable = {old + value for old in reachable for value in values if old + value <= maximum}
    return reachable


def arithmetic_admissible(value: int) -> bool:
    return valuation(value, 3) >= 14 and is_rational_eisenstein_norm(value)


def zero_internal_candidates() -> dict[str, object]:
    # Fourteen cross-color pairs cost at least 14*3=42.  Write c=42+9h.
    rows = []
    for h in range((ENERGY_LIMIT - 42) // 9 + 1):
        cross_energy = 42 + 9 * h
        determinant_value = 15**13 * (225 - cross_energy)
        if determinant_value > BENCHMARK:
            rows.append({
                "h": h,
                "cross_energy": cross_energy,
                "determinant": determinant_value,
                "v3": valuation(determinant_value, 3),
                "eisenstein_norm": is_rational_eisenstein_norm(determinant_value),
            })
    assert rows
    assert not any(arithmetic_admissible(row["determinant"]) for row in rows)
    return {
        "above_record_candidates": rows,
        "arithmetic_survivors": 0,
    }


def one_internal_edge_candidates() -> dict[str, object]:
    catalogue = inner_product_values(15)
    cross_values = tuple(entry["value"] for entry in catalogue if 0 < entry["norm"] <= 81 and entry["norm"] % 9 == 3)
    internal_values = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    cross_norms = tuple(sorted({value.norm() for value in cross_values}))
    remaining_sums = reachable_sums(cross_norms, 12, ENERGY_LIMIT)

    local_types: set[tuple[int, int]] = set()
    for u in cross_values:
        for v in cross_values:
            for g in internal_values:
                product = u * g * v.conjugate()
                numerator = 15 * (u.norm() + v.norm()) - (2 * product.a - product.b)
                local_types.add((u.norm() + v.norm(), numerator))

    candidates: set[tuple[int, int, int, int]] = set()
    for h in range((ENERGY_LIMIT - 51) // 9 + 1):
        cross_energy = 42 + 9 * h
        for pair_energy, numerator in local_types:
            singles_energy = cross_energy - pair_energy
            if singles_energy not in remaining_sums:
                continue
            # det(A)=216*15^12 and the Schur loss is the sum of the
            # twelve scalar losses plus the exact 2x2-block loss.
            schur = Fraction(15) - Fraction(singles_energy, 15) - Fraction(numerator, 216)
            determinant_value = Fraction(216 * 15**12) * schur
            if determinant_value > BENCHMARK:
                assert determinant_value.denominator == 1
                candidates.add((h, pair_energy, numerator, determinant_value.numerator))

    admissible = [row for row in sorted(candidates) if arithmetic_admissible(row[3])]
    assert candidates
    assert len(admissible) == 10
    return {
        "local_pair_types": len(local_types),
        "above_record_abstract_candidates": len(candidates),
        "distinct_above_record_determinants": len({row[3] for row in candidates}),
        "arithmetic_survivors": len(admissible),
        "surviving_determinants": sorted({row[3] for row in admissible}),
    }


def _pair_local_types(cross_values, internal_values) -> set[tuple[int, int]]:
    result = set()
    for u in cross_values:
        for v in cross_values:
            for g in internal_values:
                product = u * g * v.conjugate()
                numerator = 15 * (u.norm() + v.norm()) - (2 * product.a - product.b)
                result.add((u.norm() + v.norm(), numerator))
    return result


def two_internal_edge_candidates() -> dict[str, object]:
    """Enumerate both support graphs with two norm-9 internal edges."""

    catalogue = inner_product_values(15)
    cross_values = tuple(entry["value"] for entry in catalogue if 0 < entry["norm"] <= 81 and entry["norm"] % 9 == 3)
    internal_values = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    cross_norms = tuple(sorted({value.norm() for value in cross_values}))
    pair_types = _pair_local_types(cross_values, internal_values)

    # Two disjoint 2x2 blocks and ten scalar blocks.
    remaining_10 = reachable_sums(cross_norms, 10, ENERGY_LIMIT)
    disjoint: set[int] = set()
    for h in range((ENERGY_LIMIT - 60) // 9 + 1):
        cross_energy = 42 + 9 * h
        for energy_1, numerator_1 in pair_types:
            for energy_2, numerator_2 in pair_types:
                singles_energy = cross_energy - energy_1 - energy_2
                if singles_energy not in remaining_10:
                    continue
                schur = Fraction(15) - Fraction(singles_energy, 15) - Fraction(numerator_1 + numerator_2, 216)
                determinant_value = Fraction(216**2 * 15**10) * schur
                if determinant_value > BENCHMARK:
                    assert determinant_value.denominator == 1
                    disjoint.add(determinant_value.numerator)

    # A 3-vertex path and eleven scalar blocks.  Independent multiplication
    # of leaf rows by cube roots gauges each norm-9 edge to +3 or -3; these
    # are the two orbits of the six norm-9 catalogue values.
    path_edge_representatives = (Eisenstein(3), Eisenstein(-3))
    path_types: set[tuple[int, int]] = set()
    for g in path_edge_representatives:
        for edge in path_edge_representatives:
            for u in cross_values:
                for v in cross_values:
                    for w in cross_values:
                        # adj([[15,g,0],[g*,15,edge],[0,edge*,15]])
                        terms = (
                            216 * u.norm()
                            + 225 * v.norm()
                            + 216 * w.norm()
                            - 15 * u * g * v.conjugate()
                            - 15 * v * g.conjugate() * u.conjugate()
                            - 15 * v * edge * w.conjugate()
                            - 15 * w * edge.conjugate() * v.conjugate()
                            + u * g * edge * w.conjugate()
                            + w * edge.conjugate() * g.conjugate() * u.conjugate()
                        )
                        assert terms.b == 0
                        path_types.add((u.norm() + v.norm() + w.norm(), terms.a))

    remaining_11 = reachable_sums(cross_norms, 11, ENERGY_LIMIT)
    path: set[int] = set()
    path_determinant = 15 * (225 - 18)
    assert path_determinant == 3105
    for h in range((ENERGY_LIMIT - 60) // 9 + 1):
        cross_energy = 42 + 9 * h
        for block_energy, numerator in path_types:
            singles_energy = cross_energy - block_energy
            if singles_energy not in remaining_11:
                continue
            schur = Fraction(15) - Fraction(singles_energy, 15) - Fraction(numerator, path_determinant)
            determinant_value = Fraction(path_determinant * 15**11) * schur
            if determinant_value > BENCHMARK:
                assert determinant_value.denominator == 1
                path.add(determinant_value.numerator)

    disjoint_admissible = sorted(value for value in disjoint if arithmetic_admissible(value))
    path_admissible = sorted(value for value in path if arithmetic_admissible(value))
    return {
        "support_types": ["two_disjoint_edges", "three_vertex_path"],
        "pair_local_types": len(pair_types),
        "path_local_types": len(path_types),
        "disjoint_above_record_determinants": len(disjoint),
        "path_above_record_determinants": len(path),
        "disjoint_arithmetic_survivors": len(disjoint_admissible),
        "path_arithmetic_survivors": len(path_admissible),
        "surviving_determinants": sorted(set(disjoint_admissible + path_admissible)),
    }


def main() -> None:
    zero = zero_internal_candidates()
    one = one_internal_edge_candidates()
    two = two_internal_edge_candidates()
    print(json.dumps({
        "color_partition": [14, 1, 0],
        "leaf_internal_energy_0": zero,
        "leaf_internal_energy_9": one,
        "leaf_internal_energy_18": two,
        "theorem": "Energy 0 is excluded; energies 9 and 18 reduce to explicit finite arithmetic determinant lists.",
    }, indent=2))


if __name__ == "__main__":
    main()
