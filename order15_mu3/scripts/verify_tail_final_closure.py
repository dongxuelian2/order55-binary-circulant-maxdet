"""Final exact closure of the genuine high-energy tail Q=153,159,162,168.

This module builds on :mod:`verify_tail_structural_reduction`.  Two additional
ideas finish the two shells left there:

* a two-hub Schur--Cauchy inequality that retains the mandatory outside Gram
  entry while needing only a lower bound on ``tr(U^* A^{-1} U)``; and
* a componentwise Schur decomposition of the no-isolate 14-row block at
  Q=159, plus an exact enumeration of the 3,159 unlabeled trees on 14
  vertices in the sole connected boundary case.

For Q=153, the last all-same-color branch is closed by combining a weighted
Motzkin--Straus cap with the intertwining identity ``AH=HB``.  If there are at
most four support isolates then K6 cannot fit in the 17-unit energy budget;
if there are at least eight isolates, intertwining forces at least eight
15-eigenvalues on the opposite Gram matrix.  The three intermediate isolate
counts reduce to elementary K6 boundary supports and are handled by exact
Fischer/Schur bounds.
"""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from math import prod

import networkx as nx

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q126_boundary import capped_trace_candidates
    from order15_mu3.scripts.verify_tail_structural_reduction import (
        _sparse_blocks,
        certificate as structural_certificate,
        sqrt_upper,
    )
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, sqrt_lower
except ModuleNotFoundError:
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q126_boundary import capped_trace_candidates
    from verify_tail_structural_reduction import _sparse_blocks, certificate as structural_certificate, sqrt_upper
    from verify_trace_stability import BENCHMARK, sqrt_lower


Q153_RESIDUAL = (
    (27, 78, 48),
    (27, 87, 39),
    (36, 78, 39),
    (54, 78, 21),
)
Q159_SIZE13_RESIDUAL = ((54, 78, 27),)


def two_hub_cauchy_upper(internal: int, cross: int, outside_norm: int) -> Fraction:
    """Sharpen a sparse 13+2 Schur bound by Cauchy--Schwarz.

    Write the Schur complement as

        S = [[15-x, g-z], [conj(g-z), 15-y]],

    where ``x=u^*A^{-1}u``, ``y=v^*A^{-1}v`` and
    ``z=u^*A^{-1}v``.  Cauchy--Schwarz gives ``|z|<=sqrt(xy)``.  If
    ``s=x+y`` and ``a=|g|`` then, for fixed s,

        det(S) <= (15-s/2)^2 - max(0, a-s/2)^2.

    The right side decreases for ``0<=s<30``.  Since
    ``s >= cross/lambda_max(A)``, an exact upper spectral bracket for A gives
    a certified determinant bound.  In the five residual allocations below,
    ``s/2 < |g|`` already follows from the rational brackets, so the displayed
    expression simplifies to

        225-|g|^2-(15-|g|)s.
    """

    assert (internal, cross, outside_norm) in Q153_RESIDUAL + Q159_SIZE13_RESIDUAL
    outside_lower = sqrt_lower(Fraction(outside_norm))
    outside_upper = sqrt_upper(Fraction(outside_norm))
    values = []
    for block_determinant, _lambda_min, lambda_max in _sparse_blocks(internal):
        s_lower = Fraction(cross, 1) / lambda_max
        assert 0 < s_lower < 30
        # This certifies the first branch of the piecewise Schur--Cauchy bound.
        assert s_lower / 2 < outside_lower
        schur = (
            Fraction(225 - outside_norm)
            - (Fraction(15) - outside_upper) * s_lower
        )
        if schur > 0:
            values.append(Fraction(block_determinant) * schur)
    assert values
    return max(values)


def _component_configurations(vertices: int, units: int):
    """Connected-component abstractions with no isolated support vertices."""

    def extend(remaining_vertices: int, remaining_units: int, minimum, prefix):
        if remaining_vertices == 0 and remaining_units == 0:
            yield prefix
            return
        for size in range(2, remaining_vertices + 1):
            # A connected size-v support needs at least v-1 nonzero edges,
            # and every same-color edge costs at least one 9-energy unit.
            for energy_units in range(size - 1, remaining_units + 1):
                item = (size, energy_units)
                if item < minimum:
                    continue
                yield from extend(
                    remaining_vertices - size,
                    remaining_units - energy_units,
                    item,
                    prefix + (item,),
                )

    yield from extend(vertices, units, (2, 1), ())


def _component_lambda_upper(size: int, units: int) -> Fraction:
    # For E=A-15I, tr(E)=0 and tr(E^2)=18*units.  A zero-sum vector of
    # size v has maximum coordinate squared at most (v-1)/v of its norm.
    deviation_square = Fraction(18 * units * (size - 1), size)
    return Fraction(15) + sqrt_upper(deviation_square)


def _component_configuration_upper(configuration, cross_energy: int) -> Fraction:
    determinants = []
    lambdas = []
    correction = Fraction(0)
    for size, units in configuration:
        determinants.append(stationary_upper_dimension(size, 9 * units))
        cap = _component_lambda_upper(size, units)
        lambdas.append(cap)
        # Different-color entries are nonzero and have norm at least three.
        correction += Fraction(3 * size, 1) / cap

    # The baseline cross energy is 3*14=42.  Giving every excess unit to the
    # component with largest lambda only weakens the correction, hence is a
    # valid global lower bound without assuming a realizable distribution.
    excess = cross_energy - 42
    assert excess >= 0 and excess % 9 == 0
    if excess:
        correction += Fraction(excess, 1) / max(lambdas)

    schur = Fraction(15) - correction
    if schur <= 0:
        return Fraction(0)
    return prod(determinants, start=Fraction(1)) * schur


def _tree_matching_counts(tree: nx.Graph) -> tuple[int, ...]:
    """Exact matching counts of a rooted tree by dynamic programming."""

    root = 0
    parent = {root: -1}
    order = [root]
    for vertex in order:
        for neighbour in tree.neighbors(vertex):
            if neighbour == parent[vertex]:
                continue
            parent[neighbour] = vertex
            order.append(neighbour)

    free = {}
    total = {}
    for vertex in reversed(order):
        children = [v for v in tree.neighbors(vertex) if parent.get(v) == vertex]

        polynomial = [1]
        for child in children:
            next_polynomial = [0] * (len(polynomial) + len(total[child]) - 1)
            for i, left in enumerate(polynomial):
                for j, right in enumerate(total[child]):
                    next_polynomial[i + j] += left * right
            polynomial = next_polynomial
        free[vertex] = polynomial

        all_matchings = polynomial[:]
        for chosen in children:
            term = [1]
            for child in children:
                factor = free[child] if child == chosen else total[child]
                next_term = [0] * (len(term) + len(factor) - 1)
                for i, left in enumerate(term):
                    for j, right in enumerate(factor):
                        next_term[i + j] += left * right
                term = next_term
            if len(all_matchings) < len(term) + 1:
                all_matchings.extend([0] * (len(term) + 1 - len(all_matchings)))
            for index, value in enumerate(term):
                all_matchings[index + 1] += value
        total[vertex] = all_matchings

    return tuple(total[root])


def _tree_gram_determinant(tree: nx.Graph) -> int:
    """det(15 I + 3 A_T) from the matching polynomial of a tree."""

    order = tree.number_of_nodes()
    return sum(
        (-1) ** matching_size
        * count
        * 3 ** (2 * matching_size)
        * 15 ** (order - 2 * matching_size)
        for matching_size, count in enumerate(_tree_matching_counts(tree))
    )


@lru_cache(maxsize=1)
def _order14_tree_maximum() -> dict[str, object]:
    maximum = -1
    maximizers = []
    count = 0
    for tree in nx.nonisomorphic_trees(14):
        count += 1
        value = _tree_gram_determinant(tree)
        if value > maximum:
            maximum = value
            maximizers = [tree.copy()]
        elif value == maximum:
            maximizers.append(tree.copy())

    assert count == 3159
    assert maximum == 16802420983158456
    assert len(maximizers) == 1
    assert nx.is_isomorphic(maximizers[0], nx.path_graph(14))
    return {
        "unlabeled_trees": count,
        "maximum_determinant": maximum,
        "maximizer": "P14",
        "matching_counts": list(_tree_matching_counts(maximizers[0])),
    }


def q159_size14_certificate() -> dict[str, object]:
    rows = []
    connected_boundary = None
    for internal in (99, 108, 117):
        units = internal // 9
        cross = 159 - internal
        configurations = list(_component_configurations(14, units))
        bounds = []
        for configuration in configurations:
            upper = _component_configuration_upper(configuration, cross)
            bounds.append((configuration, upper))

        residual = [(cfg, value) for cfg, value in bounds if value >= BENCHMARK]
        if internal < 117:
            assert not residual
        else:
            assert [cfg for cfg, _value in residual] == [((14, 13),)]
            tree_data = _order14_tree_maximum()
            lambda_cap = _component_lambda_upper(14, 13)
            tree_upper = Fraction(tree_data["maximum_determinant"]) * (
                Fraction(15) - Fraction(42, 1) / lambda_cap
            )
            assert tree_upper < BENCHMARK
            connected_boundary = {
                **tree_data,
                "lambda_upper": str(lambda_cap),
                "upper_numerator": tree_upper.numerator,
                "upper_denominator": tree_upper.denominator,
                "over_record": float(tree_upper / BENCHMARK),
            }

        rows.append({
            "internal": internal,
            "cross": cross,
            "component_configurations": len(configurations),
            "generic_max_over_record": float(max(value for _cfg, value in bounds) / BENCHMARK),
            "generic_residual_configurations": [list(map(list, cfg)) for cfg, _value in residual],
        })

    assert connected_boundary is not None
    return {
        "rows": rows,
        "connected_e117_tree_boundary": connected_boundary,
        "theorem": "Every Q=159 color-(14,1) allocation is below the record.",
    }


def _q153_clique5_upper() -> Fraction:
    # If the support clique number is at most five, weighted Cauchy--Schwarz
    # plus Motzkin--Straus gives rho(E)^2 <= 2*153*(1-1/5)=1224/5.
    cap = Fraction(154, 5)
    assert Fraction(1224, 5) < (cap - 15) ** 2
    upper = max(row[3] for row in capped_trace_candidates(cap, 306))
    assert upper < BENCHMARK
    return upper


def q153_all_same_certificate() -> dict[str, object]:
    clique5_upper = _q153_clique5_upper()

    # Eight row-support isolates give eight independent rows of H in ker(B)
    # through AH=HB.  Hence the opposite Gram has at least eight eigenvalues
    # exactly 15.  The remaining seven eigenvalues carry the whole variance.
    eight_isolate_upper = Fraction(15**8) * stationary_upper_dimension(7, 153)
    assert eight_isolate_upper < BENCHMARK

    # A K6 consumes 15 of the 17 same-color energy units.  The complete
    # six-vertex principal block therefore has energy exactly 135 in every
    # intermediate-isolate boundary below.  Its fixed-trace/variance maximum
    # is 30*12^5.
    k6_upper = stationary_upper_dimension(6, 135)
    assert k6_upper == 30 * 12**5

    # m=5: K6 plus two disjoint K2 blocks and five support isolates.
    isolate5 = k6_upper * 216**2 * 15**5
    assert isolate5 < BENCHMARK

    # m=6: two remaining unit edges must cover three active vertices, so the
    # complementary 3x3 principal block contains at least one norm-9 edge.
    # Fischer then gives det <= det(K6)*(216*15)*15^6.
    isolate6 = k6_upper * (216 * 15) * 15**6
    assert isolate6 < BENCHMARK

    # m=7: if the two outside vertices are adjacent, Fischer gives the same
    # K2 bound.  Otherwise both remaining unit edges run from those vertices
    # into K6.  Their cross energy is 18; lambda_max(K6)<=30, so the 2x2
    # Schur complement has trace at most 30-18/30=147/5.
    isolate7_adjacent = k6_upper * 216 * 15**7
    isolate7_cross = k6_upper * Fraction(147, 10) ** 2 * 15**7
    assert isolate7_adjacent < BENCHMARK
    assert isolate7_cross < BENCHMARK

    return {
        "clique_at_most_five_over_record": float(clique5_upper / BENCHMARK),
        "eight_or_more_isolates_over_record": float(eight_isolate_upper / BENCHMARK),
        "k6_block_upper": k6_upper,
        "five_isolates_over_record": float(isolate5 / BENCHMARK),
        "six_isolates_over_record": float(isolate6 / BENCHMARK),
        "seven_isolates_adjacent_over_record": float(isolate7_adjacent / BENCHMARK),
        "seven_isolates_cross_over_record": float(isolate7_cross / BENCHMARK),
        "theorem": "No all-same-color Q=153 Gram matrix can beat the record.",
    }


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    prior = structural_certificate()
    assert prior["excluded_tail_shells"] == [162, 168]

    q153_improved = {
        f"{internal},{cross},{outside}": two_hub_cauchy_upper(internal, cross, outside)
        for internal, cross, outside in Q153_RESIDUAL
    }
    q159_size13_improved = {
        f"{internal},{cross},{outside}": two_hub_cauchy_upper(internal, cross, outside)
        for internal, cross, outside in Q159_SIZE13_RESIDUAL
    }
    assert all(value < BENCHMARK for value in q153_improved.values())
    assert all(value < BENCHMARK for value in q159_size13_improved.values())

    q159_size14 = q159_size14_certificate()
    q153_all_same = q153_all_same_certificate()

    return {
        "closed_tail_shells": [153, 159, 162, 168],
        "q153_size13_improved_over_record": {
            key: float(value / BENCHMARK) for key, value in q153_improved.items()
        },
        "q153_all_same": q153_all_same,
        "q159_size13_improved_over_record": {
            key: float(value / BENCHMARK) for key, value in q159_size13_improved.items()
        },
        "q159_size14": q159_size14,
        "q162": prior["q162"]["theorem"],
        "q168": prior["q168"]["theorem"],
        "new_counterexample_energy_bound": "Q <= 150",
        "theorem": (
            "No strict counterexample can lie in any genuine post-Q150 shell "
            "Q=153,159,162,168; together with trace stability at Q>=171, every "
            "strict counterexample must satisfy Q<=150."
        ),
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
