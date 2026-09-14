"""Exact Q=78 boundary enumeration for (14,1,0) with internal energy 36."""

from __future__ import annotations

import json
from itertools import product

from maxdet.mu3 import Eisenstein, determinant, inner_product_values

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_color_15_low_energy import component_multisets, connected_types
    from order15_mu3.scripts.verify_support_independence import CODE_BOUND
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_bounds import is_rational_eisenstein_norm
    from verify_color_15_low_energy import component_multisets, connected_types
    from verify_support_independence import CODE_BOUND
    from verify_trace_stability import BENCHMARK


STAR = Eisenstein(1, -1)


def gram_value(edges, labels) -> int:
    matrix = [[Eisenstein(15 if i == j else 0) for j in range(15)] for i in range(15)]
    for leaf in range(14):
        matrix[leaf][14] = STAR
        matrix[14][leaf] = STAR.conjugate()
    for (i, j), label in zip(edges, labels):
        matrix[i][j] = label
        matrix[j][i] = label.conjugate()
    value = determinant(matrix)
    assert value.b == 0 and value.a > 0
    return value.a


def enumerate_boundary() -> dict[str, object]:
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    norm27 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 27)
    types = connected_types(4, norm9)

    profiles = []
    q9_values = set()
    for indices in component_multisets(types, 4):
        components = [types[index] for index in indices]
        used = sum(component["vertices"] for component in components)
        if used > 14:
            continue
        alpha = 14 - used + sum(component["alpha"] for component in components)
        if alpha > CODE_BOUND:
            continue
        edges = []
        offset = 0
        for component in components:
            edges.extend((offset + i, offset + j) for i, j in component["edge_list"])
            offset += component["vertices"]
        values = {gram_value(edges, labels) for labels in product(norm9, repeat=4)}
        q9_values.update(values)
        profiles.append({
            "component_edges": [component["edge_list"] for component in components],
            "internal_isolates": 14 - used,
            "alpha": alpha,
            "arithmetic_survivors": sorted(value for value in values if value > BENCHMARK and is_rational_eisenstein_norm(value)),
        })

    # Weight partition 27+9.  Independence <=12 forces the two edges to be
    # disjoint; a single norm-36 edge has independence 13 and is excluded.
    mixed_edges = ((0, 1), (2, 3))
    mixed_values = {
        gram_value(mixed_edges, labels)
        for labels in product(norm27, norm9)
    } | {
        gram_value(mixed_edges, labels)
        for labels in product(norm9, norm27)
    }
    q9_survivors = sorted(value for value in q9_values if value > BENCHMARK and is_rational_eisenstein_norm(value))
    mixed_survivors = sorted(value for value in mixed_values if value > BENCHMARK and is_rational_eisenstein_norm(value))
    return {
        "color_partition": [14, 1, 0],
        "total_energy": 78,
        "internal_energy": 36,
        "cross_energy": 42,
        "four_norm9_profiles": profiles,
        "four_norm9_arithmetic_survivors": q9_survivors,
        "norm27_plus_norm9_arithmetic_survivors": mixed_survivors,
        "norm36_single_edge": "excluded by independence number 13 > 12",
    }


def main() -> None:
    result = enumerate_boundary()
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
