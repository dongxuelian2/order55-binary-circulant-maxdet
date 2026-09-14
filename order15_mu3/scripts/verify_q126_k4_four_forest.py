"""Exact phase certificate for the tight weighted K4 plus four-edge forest."""

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
    """All unlabeled K5-free graphs obtained from K4 by a four-edge forest."""

    base = nx.complete_graph(4)
    base.add_nodes_from(range(4, 8))
    optional = tuple(nx.non_edges(base))
    buckets: dict[str, list[nx.Graph]] = {}
    for added in combinations(optional, 4):
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
        graph = nx.convert_node_labels_to_integers(graph)
        edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
        result.append({
            "vertices": 8,
            "edges": 10,
            "edge_list": edges,
            "alpha": independence_number(8, edges),
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
        types,
        (1,) * 10,
        catalogue,
        required_isolates=0,
        order=8,
    )
    weighted = weighted_layer(
        types,
        (3,) + (1,) * 9,
        catalogue,
        required_isolates=0,
        order=8,
    )
    weighted4 = weighted_layer(
        types,
        (4,) + (1,) * 9,
        catalogue,
        required_isolates=0,
        order=8,
    )
    component_maxima = {
        "unit": unit["maximum_determinant"],
        "weighted": weighted["maximum_determinant"],
        "weighted4": weighted4["maximum_determinant"],
    }
    placements = {
        "heavy_outside": component_maxima["unit"] * 198 * 216 * 15**3,
        "heavy_inside": component_maxima["weighted"] * 216**2 * 15**3,
    }
    maximum = max(placements.values())
    weight4_placements = {
        "heavy_outside": component_maxima["unit"] * 189 * 15**5,
        "heavy_inside": component_maxima["weighted4"] * 216 * 15**5,
    }
    weight4_maximum = max(weight4_placements.values())
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert len(types) == 19
    assert component_maxima["unit"] == 1766431152
    assert component_maxima["weighted"] == 1628617347
    assert component_maxima["weighted4"] == 1573793631
    assert maximum == 256448601928008000
    assert weight4_maximum == 258141500324775000
    assert max(maximum, weight4_maximum) < q153
    return {
        "graph_type_count": len(types),
        "profile_counts": {
            "unit": len(unit["profiles"]),
            "weighted": len(weighted["profiles"]),
            "weighted4": len(weighted4["profiles"]),
        },
        "component_maxima": component_maxima,
        "placement_maxima": placements,
        "maximum": maximum,
        "weight4_placement_maxima": weight4_placements,
        "weight4_maximum": weight4_maximum,
        "over_q153": float(max(maximum, weight4_maximum) / q153),
        "theorem": "Both tight Q=126 weighted K4 plus four-edge-forest branches lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
