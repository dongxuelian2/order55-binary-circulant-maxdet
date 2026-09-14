"""Exact certificate excluding all-same-color total Gram energy Q=81."""

from __future__ import annotations

import json
from collections import defaultdict
from itertools import product

from pysat.solvers import Solver

from maxdet.mu3 import Eisenstein, inner_product_values

try:
    from order15_mu3.scripts.sat_q81_reduced import build as build_reduced_sat
    from order15_mu3.scripts.verify_color_15_low_energy import (
        BENCHMARK,
        connected_types,
        q9_layer,
        weighted_layer,
    )
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from sat_q81_reduced import build as build_reduced_sat
    from verify_color_15_low_energy import BENCHMARK, connected_types, q9_layer, weighted_layer


MIXED_PARTITIONS = (
    (9,),
    (7, 1, 1),
    (4, 4, 1),
    (4, 3, 1, 1),
    (4, 1, 1, 1, 1, 1),
    (3, 3, 3),
    (3, 3, 1, 1, 1),
    (3, 1, 1, 1, 1, 1, 1),
)


def collect_kernel_states(layers) -> dict[int, set[tuple[int, int, int]]]:
    """Map determinant to (isolates, active vertices, active nullity)."""

    result: dict[int, set[tuple[int, int, int]]] = defaultdict(set)
    for layer in layers:
        for profile in layer["profiles"]:
            isolates = profile["isolates"]
            if isolates == 0:
                continue
            active = 15 - isolates
            for value_text, nullities in profile["active_kernel_nullities"].items():
                value = int(value_text)
                for nullity in nullities:
                    result[value].add((isolates, active, nullity))
    return dict(result)


def rank_feasible_pairs(states_by_value) -> list[dict[str, object]]:
    """Pairs surviving both m_A v_B <= 15 d_B rank inequalities."""

    survivors = []
    for value, states in states_by_value.items():
        for left in states:
            for right in states:
                m_left, v_left, d_left = left
                m_right, v_right, d_right = right
                if m_left * v_right <= 15 * d_right and m_right * v_left <= 15 * d_left:
                    survivors.append({"determinant": value, "row_state": left, "column_state": right})
    return survivors


def triangle_balance_cases(target: Eisenstein) -> list[tuple[int, Eisenstein]]:
    """Possible outside difference count and 3-row corner contribution."""

    roots = (Eisenstein(1), Eisenstein(0, 1), Eisenstein(-1, -1))
    inside_sums = {sum((roots[index] for index in choices), Eisenstein()) for choices in product(range(3), repeat=3)}
    cases = []
    for omega_count in range(13):
        outside = omega_count * roots[1] + (12 - omega_count) * roots[2]
        if target - outside in inside_sums:
            cases.append((omega_count, target - outside))
    return cases


def main() -> None:
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)

    # Kernel-compatible isolated supports.  Q=81 is a multiple of 9 below the
    # (13,1,1) threshold 90, so the column partition is also all-same-color.
    # An isolated row supplies a third-root kernel vector on every active
    # column component.
    q9_kernel = q9_layer(connected_types(9, norm9, kernel_only=True), 9)
    mixed_types = connected_types(7, norm9, kernel_only=True)
    mixed = [weighted_layer(mixed_types, weights, catalogue, kernel_only=True) for weights in MIXED_PARTITIONS]
    states = collect_kernel_states([q9_kernel, *mixed])
    rank_survivors = rank_feasible_pairs(states)
    assert states and not rank_survivors

    # Nine q9 edges are the only way to cover all 15 vertices.  Exact
    # component enumeration leaves K3 + 6 K2 and four triangle cycle phases.
    no_isolate_layer = q9_layer(connected_types(4, norm9), 9)
    no_isolate_profiles = [profile for profile in no_isolate_layer["profiles"] if profile["isolates"] == 0]
    no_isolate_values = sorted({value for profile in no_isolate_profiles for value in profile["arithmetic_survivors"]})
    expected_values = [
        (2916 + 27 * offset) * 216**6
        for offset in (0, 1, 3, 4)
    ]
    assert no_isolate_values == expected_values
    assert not set(no_isolate_values) & set(states)

    # Triangle cycle real parts ±1/2 have no ±1 eigenvalue and contradict the
    # K2 identity h B^2=h.  For real part +1 the exact three-column balance
    # has no solution.  For real part -1 it forces the reduced block form.
    plus_cases = triangle_balance_cases(Eisenstein(3))
    minus_cases = triangle_balance_cases(Eisenstein(-3))
    assert plus_cases == []
    assert minus_cases == [(6, Eisenstein(3))]
    cnf, pool, _cell = build_reduced_sat()
    with Solver(name="glucose42", bootstrap_with=cnf) as solver:
        reduced_result = solver.solve()
    assert reduced_result is False

    print(json.dumps({
        "color_partition": [15, 0, 0],
        "energy": 81,
        "isolated_kernel_determinants": len(states),
        "isolated_kernel_states": {str(value): sorted(items) for value, items in sorted(states.items())},
        "rank_inequality_survivors": rank_survivors,
        "no_isolate_support": "K3 plus 6 K2",
        "no_isolate_arithmetic_values": no_isolate_values,
        "triangle_plus3_balance_cases": plus_cases,
        "triangle_minus3_balance_cases": [(count, [value.a, value.b]) for count, value in minus_cases],
        "reduced_sat": {"variables": pool.top, "clauses": len(cnf.clauses), "result": "unsat"},
        "theorem": "All-same-color energy Q=81 is impossible.",
    }, indent=2))


if __name__ == "__main__":
    main()
