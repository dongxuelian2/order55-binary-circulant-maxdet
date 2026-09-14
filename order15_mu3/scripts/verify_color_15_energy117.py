"""Exact all-same-color support certificate at total Gram energy Q=117."""

from __future__ import annotations

import json
from functools import lru_cache

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions
    from order15_mu3.scripts.verify_color_15_energy99 import _tree_type
    from order15_mu3.scripts.verify_color_15_low_energy import (
        connected_near_tree_graphs,
        connected_types,
        q9_layer,
        weighted_layer,
    )
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_15_energy90 import weight_partitions
    from verify_color_15_energy99 import _tree_type
    from verify_color_15_low_energy import connected_near_tree_graphs, connected_types, q9_layer, weighted_layer
    from verify_trace_stability import BENCHMARK


MIXED_PARTITIONS = tuple(partition for partition in weight_partitions(13) if partition != (1,) * 13)
NO_ISOLATE_PARTITIONS = (
    (4, 3, 1, 1, 1, 1, 1, 1),
    (3, 3, 1, 1, 1, 1, 1, 1, 1),
    (4, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    (3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1),
    (1,) * 13,
)


def balanced_unit_support_edge_divisibility() -> dict[str, object]:
    """Certificate for the pure-unit isolated-support obstruction.

    After switching a full-support third-root kernel vector to all ones, the
    incident unit gains at every vertex sum to zero.  A sum of third roots is
    zero exactly when the three root multiplicities agree, so every degree is
    divisible by three.  The degree sum is twice the edge count; hence the
    edge count itself is divisible by three.
    """

    assert 13 % 3 == 1
    return {
        "vertex_degree_modulus": 3,
        "edge_count_modulus": 3,
        "q117_unit_edges": 13,
        "q117_unit_edges_modulus": 1,
        "possible": False,
    }


def _q13_no_isolate_types(norm9):
    # If a connected component with e>8 occurs in a 13-edge support covering
    # all 15 vertices, the other 13-e edges cover at most twice that many
    # vertices.  The following near-tree pairs are exactly the possibilities.
    pairs = (
        (8, 9), (9, 9), (9, 10), (10, 9), (10, 10),
        (11, 10), (11, 11), (12, 11), (13, 12),
    )
    types = connected_types(8, norm9)
    for vertices, edges in pairs:
        types += [_tree_type(graph, norm9) for graph in connected_near_tree_graphs(vertices, edges)]
    return types


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    assert len(MIXED_PARTITIONS) == 18
    assert set(NO_ISOLATE_PARTITIONS) == {
        partition for partition in weight_partitions(13) if len(partition) >= 8
    }
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)

    # If one Gram support has an isolate, AH=HB forces the opposite active
    # support to possess a full-support mu3 kernel vector.  Exact balance
    # enumeration leaves no above-record arithmetic value for any mixed
    # weight partition.  In the pure thirteen-unit case the degree-divisibility
    # lemma eliminates the opposite support without graph enumeration.
    kernel_catalogues = {}
    isolated_mixed = {}
    for weights in MIXED_PARTITIONS:
        edge_count = len(weights)
        if edge_count not in kernel_catalogues:
            kernel_catalogues[edge_count] = connected_types(edge_count, norm9, kernel_only=True)
        layer = weighted_layer(
            kernel_catalogues[edge_count], weights, catalogue, mu3_kernel_only=True
        )
        assert not layer["intertwining_arithmetic_survivors"]
        isolated_mixed[str(weights)] = {
            "profiles": len(layer["profiles"]),
            "kernel_arithmetic_survivors": layer["intertwining_arithmetic_survivors"],
        }
    pure_isolated = balanced_unit_support_edge_divisibility()
    assert not pure_isolated["possible"]

    # With no isolates, fewer than eight support edges cannot cover 15
    # vertices.  The unique eight-edge support is P3+6K2; distributing the
    # norm-36 and norm-27 edges gives these four exact determinants.
    eight_edge_values = sorted((
        2430 * 216**6,
        2700 * 198 * 216**5,
        2835 * 189 * 216**5,
        3105 * 189 * 198 * 216**4,
    ))
    assert max(eight_edge_values) < BENCHMARK

    nine_edge = weighted_layer(
        connected_types(4, norm9), NO_ISOLATE_PARTITIONS[1], catalogue, required_isolates=0
    )
    ten_edge = weighted_layer(
        connected_types(6, norm9), NO_ISOLATE_PARTITIONS[2], catalogue, required_isolates=0
    )
    eleven_types = connected_types(7, norm9) + [
        _tree_type(graph, norm9) for graph in nx.nonisomorphic_trees(9)
    ]
    eleven_edge = weighted_layer(
        eleven_types, NO_ISOLATE_PARTITIONS[3], catalogue, required_isolates=0
    )
    assert not nine_edge["arithmetic_survivors"]
    assert not ten_edge["arithmetic_survivors"]
    assert not eleven_edge["arithmetic_survivors"]

    q13_types = _q13_no_isolate_types(norm9)
    q13_no_isolate = q9_layer(q13_types, 13, required_isolates=0)
    assert len(q13_types) == 6415
    assert len(q13_no_isolate["profiles"]) == 10329
    assert not q13_no_isolate["arithmetic_survivors"]

    return {
        "energy": 117,
        "mixed_weight_partitions": MIXED_PARTITIONS,
        "isolated_mixed_counts": isolated_mixed,
        "isolated_pure_unit_obstruction": pure_isolated,
        "no_isolate_eight_edge_values": eight_edge_values,
        "no_isolate_nine_edge_arithmetic_values": nine_edge["arithmetic_survivors"],
        "no_isolate_ten_edge_arithmetic_values": ten_edge["arithmetic_survivors"],
        "no_isolate_eleven_edge_arithmetic_values": eleven_edge["arithmetic_survivors"],
        "no_isolate_q13_component_types": len(q13_types),
        "no_isolate_q13_profiles": len(q13_no_isolate["profiles"]),
        "no_isolate_q13_arithmetic_values": q13_no_isolate["arithmetic_survivors"],
        "theorem": "No all-same-color Q=117 Gram matrix can beat the benchmark.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
