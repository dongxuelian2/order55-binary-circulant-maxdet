"""Exact phase certificate for the three tight active-order-five triangle rows."""

from __future__ import annotations

import json
from functools import lru_cache

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_low_energy import connected_near_tree_graphs, weighted_layer
    from order15_mu3.scripts.verify_support_independence import independence_number
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_color_15_low_energy import connected_near_tree_graphs, weighted_layer
    from verify_support_independence import independence_number
    from verify_trace_stability import stationary_upper


def _types(edge_count: int) -> tuple[dict[str, object], ...]:
    result = []
    for graph in connected_near_tree_graphs(5, edge_count):
        if not any(len(cycle) == 3 for cycle in nx.cycle_basis(graph)):
            continue
        if max(map(len, nx.find_cliques(graph))) >= 4:
            continue
        edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
        result.append({
            "vertices": 5,
            "edges": edge_count,
            "edge_list": edges,
            "alpha": independence_number(5, edges),
            "has_leaf": any(graph.degree(vertex) == 1 for vertex in graph),
            "determinants": set(),
            "states": set(),
        })
    return tuple(result)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    catalogue = inner_product_values(15)
    specifications = {
        "edges10_component5": (5, (3, 3)),
        "edges11_component6": (6, (4,)),
        "edges12_component7": (7, (3,)),
    }
    graph_counts = {}
    allocation_maxima = {}
    maxima = {}
    for name, (component_edges, heavy_weights) in specifications.items():
        types = _types(component_edges)
        graph_counts[name] = len(types)
        rows = {}
        for heavy_inside in range(len(heavy_weights) + 1):
            weights = heavy_weights[:heavy_inside] + (1,) * (component_edges - heavy_inside)
            layer = weighted_layer(
                types, weights, catalogue, required_isolates=0, order=5
            )
            outside_heavy = len(heavy_weights) - heavy_inside
            heavy_k2 = 189 if heavy_weights == (4,) else 198
            outside = (heavy_k2,) * outside_heavy + (216,) * (5 - outside_heavy)
            rows[str(heavy_inside)] = layer["maximum_determinant"] * __import__("math").prod(outside)
        allocation_maxima[name] = rows
        maxima[name] = max(rows.values())
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert graph_counts == {
        "edges10_component5": 3,
        "edges11_component6": 4,
        "edges12_component7": 3,
    }
    assert maxima == {
        "edges10_component5": 250675363046817792,
        "edges11_component6": 254331521486880768,
        "edges12_component7": 255359816048148480,
    }
    assert all(value < q153 for value in maxima.values())
    return {
        "graph_type_counts": graph_counts,
        "allocation_maxima": allocation_maxima,
        "maxima": maxima,
        "maximum": max(maxima.values()),
        "over_q153": float(max(maxima.values()) / q153),
        "theorem": "All three tight active-order-five Q=126 triangle rows lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
