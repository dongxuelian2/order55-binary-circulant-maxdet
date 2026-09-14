"""Exact phase certificate for weighted Q=126 K4 plus a three-edge forest."""

from __future__ import annotations

import json
from functools import lru_cache
from itertools import combinations

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_low_energy import weighted_layer
    from order15_mu3.scripts.verify_support_independence import independence_number
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_color_15_low_energy import weighted_layer
    from verify_support_independence import independence_number
    from verify_trace_stability import stationary_upper


def _support_types() -> tuple[dict[str, object], ...]:
    base = nx.complete_graph(4)
    base.add_nodes_from(range(4, 7))
    buckets: dict[str, list[nx.Graph]] = {}
    for added in combinations(tuple(nx.non_edges(base)), 3):
        graph = base.copy()
        graph.add_edges_from(added)
        if not nx.is_connected(graph) or max(map(len, nx.find_cliques(graph))) >= 5:
            continue
        key = nx.weisfeiler_lehman_graph_hash(graph)
        bucket = buckets.setdefault(key, [])
        if not any(nx.is_isomorphic(graph, prior) for prior in bucket):
            bucket.append(graph)
    result = []
    for graph in (item for bucket in buckets.values() for item in bucket):
        edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
        result.append({
            "vertices": 7,
            "edges": 9,
            "edge_list": edges,
            "alpha": independence_number(7, edges),
            "has_leaf": any(graph.degree(vertex) == 1 for vertex in graph),
            "determinants": set(),
            "states": set(),
        })
    return tuple(result)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    types = _support_types()
    catalogue = inner_product_values(15)
    unit = weighted_layer(
        types, (1,) * 9, catalogue, required_isolates=0, order=7
    )
    weighted = weighted_layer(
        types, (3,) + (1,) * 8, catalogue, required_isolates=0, order=7
    )
    component_maxima = {
        "unit": unit["maximum_determinant"],
        "weighted": weighted["maximum_determinant"],
    }
    placements = {
        "heavy_outside": component_maxima["unit"] * 198 * 216**2 * 15**2,
        "heavy_inside": component_maxima["weighted"] * 216**3 * 15**2,
    }
    maximum = max(placements.values())
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert len(types) == 7
    assert component_maxima == {"unit": 122891904, "weighted": 113304096}
    assert maximum == 256914952884633600
    assert maximum < q153
    return {
        "graph_type_count": len(types),
        "profile_counts": {
            "unit": len(unit["profiles"]),
            "weighted": len(weighted["profiles"]),
        },
        "component_maxima": component_maxima,
        "placement_maxima": placements,
        "maximum": maximum,
        "over_q153": float(maximum / q153),
        "theorem": "The tight weighted Q=126 K4 plus three-edge forest lies below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
