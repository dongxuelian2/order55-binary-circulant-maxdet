"""Exact all-same-color support certificate at total Gram energy Q=108."""

from __future__ import annotations

import json
from functools import lru_cache

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy81 import collect_kernel_states, rank_feasible_pairs
    from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions
    from order15_mu3.scripts.verify_color_15_energy99 import _tree_type
    from order15_mu3.scripts.verify_color_15_low_energy import (
        connected_types,
        connected_unicyclic_graphs,
        q9_layer,
        weighted_layer,
    )
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_15_energy81 import collect_kernel_states, rank_feasible_pairs
    from verify_color_15_energy90 import weight_partitions
    from verify_color_15_energy99 import _tree_type
    from verify_color_15_low_energy import connected_types, connected_unicyclic_graphs, q9_layer, weighted_layer
    from verify_trace_stability import BENCHMARK


MIXED_PARTITIONS = tuple(partition for partition in weight_partitions(12) if partition != (1,) * 12)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    assert len(MIXED_PARTITIONS) == 15
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
    q12_isolated = q9_layer(connected_types(12, norm9, kernel_only=True), 12)
    isolated_layers.append(q12_isolated)
    states = collect_kernel_states(isolated_layers)
    rank_survivors = rank_feasible_pairs(states)
    assert states and not rank_survivors

    # With no isolates, only the 8-, 9-, 10-, and 12-edge partitions remain.
    eight_edge_values = sorted((2565 * 216**6, 2835 * 198 * 216**5, 3105 * 198**2 * 216**4))
    assert max(eight_edge_values) < BENCHMARK
    nine_edge = weighted_layer(
        connected_types(4, norm9), (4,) + (1,) * 8, catalogue, required_isolates=0
    )
    ten_edge = weighted_layer(
        connected_types(6, norm9), (3,) + (1,) * 9, catalogue, required_isolates=0
    )
    assert not nine_edge["arithmetic_survivors"]
    assert not ten_edge["arithmetic_survivors"]

    # A component with more than eight of the twelve edges can cover enough
    # vertices only as a 9-vertex unicyclic graph, a 10-vertex tree, or an
    # 11-vertex tree.  These are the exact types missing from the atlas-based
    # catalogue through eight edges.
    q12_types = connected_types(8, norm9)
    q12_types += [_tree_type(graph, norm9) for graph in connected_unicyclic_graphs(9)]
    for vertices in (10, 11):
        q12_types += [_tree_type(graph, norm9) for graph in nx.nonisomorphic_trees(vertices)]
    q12_no_isolate = q9_layer(q12_types, 12, required_isolates=0)
    assert not q12_no_isolate["arithmetic_survivors"]
    return {
        "energy": 108,
        "mixed_weight_partitions": MIXED_PARTITIONS,
        "isolated_mixed_counts": mixed_counts,
        "isolated_q12_profiles": len(q12_isolated["profiles"]),
        "isolated_kernel_states": {str(value): sorted(items) for value, items in sorted(states.items())},
        "isolated_rank_pair_survivors": rank_survivors,
        "no_isolate_eight_edge_values": eight_edge_values,
        "no_isolate_nine_edge_arithmetic_values": nine_edge["arithmetic_survivors"],
        "no_isolate_ten_edge_arithmetic_values": ten_edge["arithmetic_survivors"],
        "no_isolate_q12_profiles": len(q12_no_isolate["profiles"]),
        "no_isolate_q12_arithmetic_values": q12_no_isolate["arithmetic_survivors"],
        "theorem": "No all-same-color Q=108 Gram matrix can beat the benchmark.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
