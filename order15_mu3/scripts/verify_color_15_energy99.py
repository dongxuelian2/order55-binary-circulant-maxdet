"""Exact all-same-color support certificate at total Gram energy Q=99."""

from __future__ import annotations

import json
from functools import lru_cache

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy81 import collect_kernel_states, rank_feasible_pairs
    from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions
    from order15_mu3.scripts.verify_color_15_low_energy import (
        component_states,
        connected_types,
        q9_layer,
        weighted_layer,
    )
    from order15_mu3.scripts.verify_support_independence import independence_number
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_15_energy81 import collect_kernel_states, rank_feasible_pairs
    from verify_color_15_energy90 import weight_partitions
    from verify_color_15_low_energy import component_states, connected_types, q9_layer, weighted_layer
    from verify_support_independence import independence_number
    from verify_trace_stability import BENCHMARK


MIXED_PARTITIONS = tuple(partition for partition in weight_partitions(11) if partition != (1,) * 11)
Q99_ALL_SAME_UPPER = 282438239494864896


def _tree_type(graph: nx.Graph, norm9) -> dict[str, object]:
    graph = nx.convert_node_labels_to_integers(graph)
    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    states = component_states(graph, norm9)
    return {
        "vertices": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "edge_list": edges,
        "alpha": independence_number(graph.number_of_nodes(), edges),
        "has_leaf": True,
        "determinants": {state[0] for state in states},
        "states": states,
    }


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    assert len(MIXED_PARTITIONS) == 12
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)

    isolated_layers = []
    mixed_counts = {}
    for weights in MIXED_PARTITIONS:
        layer = weighted_layer(
            connected_types(len(weights), norm9, kernel_only=True),
            weights,
            catalogue,
            mu3_kernel_only=True,
        )
        isolated_layers.append(layer)
        mixed_counts[str(weights)] = {
            "profiles": len(layer["profiles"]),
            "kernel_arithmetic_survivors": layer["intertwining_arithmetic_survivors"],
        }
    q9_isolated = q9_layer(connected_types(11, norm9, kernel_only=True), 11)
    isolated_layers.append(q9_isolated)
    states = collect_kernel_states(isolated_layers)
    rank_survivors = rank_feasible_pairs(states)
    assert states and not rank_survivors

    # No isolates require at least eight support edges.  Hence only the
    # partitions 36+7*9, 27+8*9, and 11*9 remain.
    eight_edge_values = sorted((2700 * 216**6, 3105 * 189 * 216**5))
    assert max(eight_edge_values) < BENCHMARK
    nine_edge = weighted_layer(
        connected_types(4, norm9),
        (3,) + (1,) * 8,
        catalogue,
        required_isolates=0,
    )
    assert nine_edge["arithmetic_survivors"] == [Q99_ALL_SAME_UPPER]

    # For eleven unit edges, a component with eight edges can cover enough
    # vertices only when it is a nine-vertex tree; all other components have
    # at most seven edges.  Add exactly those 47 missing atlas types.
    q11_types = connected_types(7, norm9) + [
        _tree_type(graph, norm9) for graph in nx.nonisomorphic_trees(9)
    ]
    q11_no_isolate = q9_layer(q11_types, 11, required_isolates=0)
    expected_q11 = [
        278025142002757632,
        278527114981945152,
        278569241062723125,
        279539020106846208,
        279924630567321600,
        280017198236160000,
        280381650372329472,
        281036368270070784,
        281078873865090000,
        282057389657358336,
        282438239494864896,
    ]
    assert q11_no_isolate["arithmetic_survivors"] == expected_q11
    return {
        "energy": 99,
        "mixed_weight_partitions": MIXED_PARTITIONS,
        "isolated_mixed_counts": mixed_counts,
        "isolated_q9_profiles": len(q9_isolated["profiles"]),
        "isolated_kernel_states": {str(value): sorted(items) for value, items in sorted(states.items())},
        "isolated_rank_pair_survivors": rank_survivors,
        "no_isolate_eight_edge_values": eight_edge_values,
        "no_isolate_nine_edge_arithmetic_values": nine_edge["arithmetic_survivors"],
        "no_isolate_q11_profiles": len(q11_no_isolate["profiles"]),
        "no_isolate_q11_arithmetic_values": expected_q11,
        "all_same_color_upper": Q99_ALL_SAME_UPPER,
        "theorem": "Every all-same-color Q=99 Gram matrix has determinant at most 282438239494864896.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
