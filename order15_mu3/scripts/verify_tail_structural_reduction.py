"""Exact structural reduction of the genuine Q=153,159,162,168 tail.

This strengthens the shellwise trace envelope by combining:

* color-partition congruence,
* exact abstract-Gram enumeration for sparse 13-row internal blocks,
* exact rational spectral brackets certified by Sylvester's criterion,
* two-hub Schur residuals retaining the mandatory outside Gram entry,
* exact internal determinant maxima through energy 72, and
* weighted Motzkin--Straus spectral caps for the remaining dense slices.

The main new conclusion is that Q=162 and Q=168 are impossible for a strict
counterexample.  The post-Q150 tail therefore reduces to Q in {153,159}.
"""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from itertools import permutations, product
from math import ceil, comb, isqrt

import networkx as nx

from maxdet.mu3 import Eisenstein, UNITS, determinant, inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions
    from order15_mu3.scripts.verify_color_15_low_energy import (
        component_multisets,
        connected_types,
        weighted_layer,
    )
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q126_boundary import capped_trace_candidates
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, sqrt_lower
except ModuleNotFoundError:
    from verify_color_15_energy90 import weight_partitions
    from verify_color_15_low_energy import component_multisets, connected_types, weighted_layer
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q126_boundary import capped_trace_candidates
    from verify_trace_stability import BENCHMARK, sqrt_lower


SPECTRAL_GRID = 20  # exact 1/20-spaced spectral brackets
SPARSE_INTERNAL_ENERGIES = (9, 18, 27, 36, 45, 54)
EXPECTED_EXACT_INTERNAL_MAXIMA = {
    36: 1674039150000000,
    45: 1619234296875000,
    54: 1594323000000000,
    63: 1530550080000000,
    72: 1474748775000000,
}


def sqrt_upper(value: Fraction, scale: int = 10**6) -> Fraction:
    """Return a certified rational upper bound for a positive square root."""

    if value <= 0:
        return Fraction(0)
    quotient = value.numerator * scale * scale // value.denominator
    numerator = isqrt(quotient)
    result = Fraction(numerator, scale)
    if result * result < value:
        result += Fraction(1, scale)
    assert result * result >= value
    return result


def _catalogue():
    return inner_product_values(15)


def _same_color_values_by_unit() -> dict[int, tuple[Eisenstein, ...]]:
    catalogue = _catalogue()
    return {
        unit: tuple(entry["value"] for entry in catalogue if entry["norm"] == 9 * unit)
        for unit in (1, 3, 4, 7, 9)
    }


def _orbit_representatives(values: tuple[Eisenstein, ...]) -> tuple[Eisenstein, ...]:
    remaining = set(values)
    representatives = []
    while remaining:
        representative = min(remaining, key=lambda value: (value.a, value.b))
        orbit = {representative * unit for unit in UNITS}
        representatives.append(representative)
        remaining.difference_update(orbit)
    return tuple(representatives)


def _positive_definite_shift(
    block: tuple[tuple[Eisenstein, ...], ...], shift: Fraction, lower: bool
) -> bool:
    """Certify A-shift*I>0 (lower=True) or shift*I-A>0 exactly."""

    numerator = shift.numerator
    denominator = shift.denominator
    size = len(block)
    for leading in range(1, size + 1):
        matrix = []
        for i in range(leading):
            row = []
            for j in range(leading):
                if lower:
                    value = denominator * block[i][j]
                    if i == j:
                        value -= numerator
                else:
                    value = -denominator * block[i][j]
                    if i == j:
                        value += numerator
                row.append(value)
            matrix.append(row)
        minor = determinant(matrix)
        assert minor.b == 0
        if minor.a <= 0:
            return False
    return True


def _spectral_bracket(block: tuple[tuple[Eisenstein, ...], ...]) -> tuple[Fraction, Fraction]:
    """Return exact lower/upper eigenvalue brackets on a 1/20 grid."""

    # A positive-definite principal Gram block is mandatory for a strict
    # counterexample.  Abstract phase states failing this test can be dropped.
    if not _positive_definite_shift(block, Fraction(0), True):
        return Fraction(0), Fraction(0)

    denominator = SPECTRAL_GRID

    # Largest grid point L with A-LI positive definite.
    lo, hi = 0, 15 * denominator
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if _positive_definite_shift(block, Fraction(mid, denominator), True):
            lo = mid
        else:
            hi = mid
    lower = Fraction(lo, denominator)

    # Smallest grid point U with UI-A positive definite.
    lo, hi = 15 * denominator, 45 * denominator
    assert _positive_definite_shift(block, Fraction(hi, denominator), False)
    while lo + 1 < hi:
        mid = (lo + hi) // 2
        if _positive_definite_shift(block, Fraction(mid, denominator), False):
            hi = mid
        else:
            lo = mid
    upper = Fraction(hi, denominator)

    assert lower > 0 and lower < upper
    return lower, upper


@lru_cache(maxsize=None)
def _sparse_blocks(internal_energy: int) -> tuple[tuple[int, Fraction, Fraction], ...]:
    """Enumerate abstract sparse 13-row internal Gram blocks exactly."""

    assert internal_energy in SPARSE_INTERNAL_ENERGIES
    total_units = internal_energy // 9
    catalogue = _catalogue()
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    values_by_unit = _same_color_values_by_unit()
    result: set[tuple[int, Fraction, Fraction]] = set()

    for weights in weight_partitions(total_units):
        edge_count = len(weights)
        types = connected_types(edge_count, norm9)
        representatives = {
            unit: _orbit_representatives(values_by_unit[unit]) for unit in set(weights)
        }
        weight_orders = tuple(sorted(set(permutations(weights))))

        for indices in component_multisets(types, edge_count):
            components = [types[index] for index in indices]
            used_vertices = sum(component["vertices"] for component in components)
            if used_vertices > 13:
                continue
            isolates = 13 - used_vertices

            edges: list[tuple[int, int]] = []
            tree_edges: set[tuple[int, int]] = set()
            offset = 0
            for component in components:
                component_edges = tuple(component["edge_list"])
                edges.extend((offset + i, offset + j) for i, j in component_edges)
                graph = nx.Graph()
                graph.add_nodes_from(range(component["vertices"]))
                graph.add_edges_from(component_edges)
                for i, j in nx.minimum_spanning_tree(graph).edges():
                    tree_edges.add((offset + min(i, j), offset + max(i, j)))
                offset += component["vertices"]

            for ordered_weights in weight_orders:
                options = tuple(
                    representatives[unit] if edge in tree_edges else values_by_unit[unit]
                    for edge, unit in zip(edges, ordered_weights)
                )
                for labels in product(*options):
                    active = [
                        [Eisenstein(15 if i == j else 0) for j in range(used_vertices)]
                        for i in range(used_vertices)
                    ]
                    for (i, j), label in zip(edges, labels):
                        active[i][j] = label
                        active[j][i] = label.conjugate()
                    block = tuple(tuple(row) for row in active)
                    lower, upper = _spectral_bracket(block)
                    if lower == 0:
                        continue
                    active_determinant = determinant(block)
                    assert active_determinant.b == 0 and active_determinant.a > 0
                    full_determinant = active_determinant.a * 15**isolates
                    result.add((full_determinant, lower, upper))

    assert result
    return tuple(sorted(result))


def _cross_column_product(cross_energy: int) -> int:
    """Largest E1*E2 for E1+E2=c and Ei=39 mod 9, Ei>=39."""

    pairs = [
        (left, cross_energy - left)
        for left in range(39, cross_energy - 38, 9)
        if cross_energy - left >= 39 and (cross_energy - left - 39) % 9 == 0
    ]
    assert pairs
    return max(left * right for left, right in pairs)


def _sparse_allocation_upper(internal: int, cross: int, outside_norm: int) -> Fraction:
    """Exact two-hub Schur upper bound for one sparse size-13 allocation."""

    product_upper = sqrt_upper(Fraction(_cross_column_product(cross)))
    outside_lower = sqrt_lower(Fraction(outside_norm))
    values = []
    for block_determinant, lambda_min, lambda_max in _sparse_blocks(internal):
        trace_upper = Fraction(30) - Fraction(cross, 1) / lambda_max
        correction_upper = product_upper / lambda_min
        residual_lower = max(Fraction(0), outside_lower - correction_upper)
        schur_upper = (trace_upper / 2) ** 2 - residual_lower**2
        if schur_upper > 0:
            values.append(Fraction(block_determinant) * schur_upper)
    return max(values, default=Fraction(0))


@lru_cache(maxsize=None)
def _exact_internal_maximum(internal_energy: int) -> int:
    """Exact abstract 13-row internal determinant maximum for e=63 or 72."""

    assert internal_energy in (63, 72)
    total_units = internal_energy // 9
    catalogue = _catalogue()
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    types = connected_types(total_units, norm9)
    maximum = max(
        weighted_layer(types, weights, catalogue, order=13)["maximum_determinant"]
        for weights in weight_partitions(total_units)
    )
    assert maximum == EXPECTED_EXACT_INTERNAL_MAXIMA[internal_energy]
    return maximum


def _size13_dense_upper(internal: int, cross: int) -> Fraction:
    if internal == 63:
        return Fraction(_exact_internal_maximum(63)) * (
            Fraction(15) - Fraction(cross, 2 * 25)
        ) ** 2
    if internal == 72:
        return Fraction(_exact_internal_maximum(72)) * (
            Fraction(15) - Fraction(cross, 2) / Fraction(127, 5)
        ) ** 2
    if internal == 81:
        # Nine support-weight units imply at most nine support edges, hence no K5.
        # Weighted Motzkin--Straus gives rho^2 <= 2*81*(1-1/4)=243/2.
        cap = Fraction(261, 10)
        assert Fraction(243, 2) < (cap - 15) ** 2
        return stationary_upper_dimension(13, 81) * (
            Fraction(15) - Fraction(cross, 2) / cap
        ) ** 2
    if internal == 90:
        # Ten support-weight units imply clique number at most five.
        # rho^2 <= 2*90*(1-1/5)=144, so lambda_max <= 27.
        cap = Fraction(27)
        return stationary_upper_dimension(13, 90) * (
            Fraction(15) - Fraction(cross, 2) / cap
        ) ** 2
    raise AssertionError(internal)


def _outside_norms(q: int) -> tuple[int, ...]:
    catalogue = _catalogue()
    residue = 3 if q % 9 == 0 else 0
    return tuple(sorted({
        int(entry["norm"])
        for entry in catalogue
        if int(entry["norm"]) <= 81 and int(entry["norm"]) % 9 == residue
    }))


def _size13_rows(q: int) -> list[dict[str, object]]:
    rows = []
    for internal in range(9, 91, 9):
        for outside in _outside_norms(q):
            cross = q - internal - outside
            if cross < 78 or cross % 9 != 6:
                continue
            if internal <= 54:
                upper = _sparse_allocation_upper(internal, cross, outside)
                method = "exact sparse block + rational spectral Schur residual"
            else:
                upper = _size13_dense_upper(internal, cross)
                method = "exact/internal stationary block + support spectral cap"
            rows.append({
                "internal": internal,
                "cross": cross,
                "outside_norm": outside,
                "upper_numerator": upper.numerator,
                "upper_denominator": upper.denominator,
                "over_record": float(upper / BENCHMARK),
                "excluded": upper < BENCHMARK,
                "method": method,
            })
    return rows


def _all_same_upper(q: int) -> Fraction:
    assert q in (153, 162)
    # q/9 <= 18 support-weight units; K7 already costs 21 unit edges, so
    # clique number <=6.  Weighted Motzkin--Straus gives rho^2<=2Q*5/6.
    rho_square = Fraction(2 * q * 5, 6)
    cap = Fraction(31) if q == 153 else Fraction(63, 2)
    assert rho_square < (cap - 15) ** 2
    return max(row[3] for row in capped_trace_candidates(cap, 2 * q))


def _size14_rows_q159() -> list[dict[str, object]]:
    """Reduce Q=159, type (14,1), using the prior no-isolate lemma."""

    rows = []
    for internal in range(18, 118, 9):
        cross = 159 - internal
        if cross < 42 or cross % 9 != 6:
            continue
        units = internal // 9
        # An internal isolate is already excluded by the equal-color-sign
        # orbit lemma.  With <=13 support-weight units, a K5 plus enough edges
        # to cover the other nine vertices is impossible; hence omega<=4.
        assert comb(5, 2) + ceil((14 - 5) / 2) > units
        cap = Fraction(15) + sqrt_upper(Fraction(3 * internal, 2))
        upper = stationary_upper_dimension(14, internal) * (
            Fraction(15) - Fraction(cross, 1) / cap
        )
        rows.append({
            "internal": internal,
            "cross": cross,
            "lambda_cap": str(cap),
            "upper_numerator": upper.numerator,
            "upper_denominator": upper.denominator,
            "over_record": float(upper / BENCHMARK),
            "excluded": upper < BENCHMARK,
        })
    return rows


def _size14_rows_q168() -> list[dict[str, object]]:
    rows = []
    for internal in range(18, 127, 9):
        cross = 168 - internal
        if cross < 42 or cross % 9 != 6:
            continue
        # Generic trace-energy control of the largest eigenvalue of the
        # fourteen-dimensional internal block.
        cap = Fraction(15) + sqrt_upper(Fraction(2 * internal * 13, 14))
        upper = stationary_upper_dimension(14, internal) * (
            Fraction(15) - Fraction(cross, 1) / cap
        )
        rows.append({
            "internal": internal,
            "cross": cross,
            "lambda_cap": str(cap),
            "over_record": float(upper / BENCHMARK),
            "excluded": upper < BENCHMARK,
        })
    return rows


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    # Replay the sparse exact maxima used by the tail reduction.
    sparse_maxima = {
        energy: max(value[0] for value in _sparse_blocks(energy))
        for energy in (36, 45, 54)
    }
    assert sparse_maxima == {
        energy: EXPECTED_EXACT_INTERNAL_MAXIMA[energy] for energy in (36, 45, 54)
    }
    exact_maxima = {
        **sparse_maxima,
        63: _exact_internal_maximum(63),
        72: _exact_internal_maximum(72),
    }

    q153_rows = _size13_rows(153)
    q159_rows = _size13_rows(159)
    q162_rows = _size13_rows(162)
    q168_rows = _size13_rows(168)

    q153_all_same = _all_same_upper(153)
    q162_all_same = _all_same_upper(162)
    q159_size14 = _size14_rows_q159()
    q168_size14 = _size14_rows_q168()

    assert q153_all_same > BENCHMARK
    assert q162_all_same < BENCHMARK
    assert all(row["excluded"] for row in q162_rows)
    assert all(row["excluded"] for row in q168_rows)
    assert all(row["excluded"] for row in q168_size14)

    q153_residual = [row for row in q153_rows if not row["excluded"]]
    q159_size13_residual = [row for row in q159_rows if not row["excluded"]]
    q159_size14_residual = [row for row in q159_size14 if not row["excluded"]]

    assert {
        (row["internal"], row["cross"], row["outside_norm"])
        for row in q153_residual
    } == {
        (27, 78, 48),
        (27, 87, 39),
        (36, 78, 39),
        (54, 78, 21),
    }
    assert {
        (row["internal"], row["cross"], row["outside_norm"])
        for row in q159_size13_residual
    } == {(54, 78, 27)}
    assert {row["internal"] for row in q159_size14_residual} == {99, 108, 117}

    return {
        "exact_size13_internal_maxima": exact_maxima,
        "excluded_tail_shells": [162, 168],
        "remaining_tail_shells": [153, 159],
        "q153": {
            "all_same_over_record": float(q153_all_same / BENCHMARK),
            "size13_rows": q153_rows,
            "residual_size13_allocations": q153_residual,
            "residual_branch_count": 1 + len(q153_residual),
        },
        "q159": {
            "size13_rows": q159_rows,
            "size14_rows": q159_size14,
            "residual_size13_allocations": q159_size13_residual,
            "residual_size14_allocations": q159_size14_residual,
            "residual_branch_count": len(q159_size13_residual) + len(q159_size14_residual),
        },
        "q162": {
            "all_same_over_record": float(q162_all_same / BENCHMARK),
            "size13_rows": q162_rows,
            "theorem": "Q=162 cannot contain a strict counterexample.",
        },
        "q168": {
            "size13_rows": q168_rows,
            "size14_rows": q168_size14,
            "theorem": "Q=168 cannot contain a strict counterexample.",
        },
        "theorem": (
            "After the exact Q<=150 branch certificates, a strict counterexample "
            "in the high-energy tail can only have Q=153 or Q=159."
        ),
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
