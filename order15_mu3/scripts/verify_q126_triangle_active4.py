"""Exact phase certificate for the tight active-order-four triangle rows."""

from __future__ import annotations

import json
from functools import lru_cache
from math import prod

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
    for graph in connected_near_tree_graphs(4, edge_count):
        if not any(len(cycle) == 3 for cycle in nx.cycle_basis(graph)):
            continue
        if max(map(len, nx.find_cliques(graph))) >= 4:
            continue
        edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
        result.append({
            "vertices": 4,
            "edges": edge_count,
            "edge_list": edges,
            "alpha": independence_number(4, edges),
            "has_leaf": any(graph.degree(vertex) == 1 for vertex in graph),
            "determinants": set(),
            "states": set(),
        })
    return tuple(result)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    catalogue = inner_product_values(15)
    specifications = {
        "edges9_component4": (4, (4, 3)),
        "edges10_component5": (5, (3, 3)),
    }
    graph_counts = {}
    allocation_maxima = {}
    maxima = {}
    for name, (component_edges, heavy_weights) in specifications.items():
        types = _types(component_edges)
        graph_counts[name] = len(types)
        rows = {}
        allocations = (
            ((inside4, inside3) for inside4 in (0, 1) for inside3 in (0, 1))
            if heavy_weights == (4, 3)
            else ((0, count) for count in range(3))
        )
        for inside4, inside3 in allocations:
            inside = (4,) * inside4 + (3,) * inside3
            weights = inside + (1,) * (component_edges - len(inside))
            layer = weighted_layer(
                types, weights, catalogue, required_isolates=0, order=4
            )
            outside = (
                (189,) * (heavy_weights.count(4) - inside4)
                + (198,) * (heavy_weights.count(3) - inside3)
                + (216,) * (5 - len(heavy_weights) + inside4 + inside3)
            )
            rows[f"{inside4},{inside3}"] = layer["maximum_determinant"] * prod(outside) * 15
        allocation_maxima[name] = rows
        maxima[name] = max(rows.values())
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert graph_counts == {"edges9_component4": 1, "edges10_component5": 1}
    assert maxima == {
        "edges9_component4": 248218881594900480,
        "edges10_component5": 254788541291888640,
    }
    assert all(value < q153 for value in maxima.values())
    return {
        "graph_type_counts": graph_counts,
        "allocation_maxima": allocation_maxima,
        "maxima": maxima,
        "maximum": max(maxima.values()),
        "over_q153": float(max(maxima.values()) / q153),
        "theorem": "Both tight active-order-four Q=126 triangle rows lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
