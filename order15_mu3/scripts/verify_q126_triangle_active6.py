"""Exact phase certificate for the three tight active-order-six triangle rows."""

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
    graphs = (
        graph for graph in connected_near_tree_graphs(6, edge_count)
        if any(len(cycle) == 3 for cycle in nx.cycle_basis(graph))
        and max(map(len, nx.find_cliques(graph))) < 4
    )
    result = []
    for graph in graphs:
        edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
        result.append({
            "vertices": 6,
            "edges": edge_count,
            "edge_list": edges,
            "alpha": independence_number(6, edges),
            "has_leaf": any(graph.degree(vertex) == 1 for vertex in graph),
            "determinants": set(),
            "states": set(),
        })
    return tuple(result)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    catalogue = inner_product_values(15)
    specifications = {
        "edges10_component6": (6, (3, 3), (198, 198, 216, 216)),
        "edges11_component7": (7, (4,), (189, 216, 216, 216)),
        "edges12_component8": (8, (3,), (198, 216, 216, 216)),
    }
    graph_counts = {}
    maxima = {}
    allocation_maxima = {}
    for name, (component_edges, heavy_weights, outside_default) in specifications.items():
        types = _types(component_edges)
        graph_counts[name] = len(types)
        rows = {}
        for heavy_inside in range(len(heavy_weights) + 1):
            inside = heavy_weights[:heavy_inside] + (1,) * (component_edges - heavy_inside)
            layer = weighted_layer(
                types, inside, catalogue, required_isolates=0, order=6
            )
            outside_heavy = len(heavy_weights) - heavy_inside
            if heavy_weights == (3, 3):
                outside = (198,) * outside_heavy + (216,) * (4 - outside_heavy)
            elif heavy_weights == (4,):
                outside = ((189,) if outside_heavy else ()) + (216,) * (4 - outside_heavy)
            else:
                outside = ((198,) if outside_heavy else ()) + (216,) * (4 - outside_heavy)
            rows[str(heavy_inside)] = layer["maximum_determinant"] * __import__("math").prod(outside) * 15
        allocation_maxima[name] = rows
        maxima[name] = max(rows.values())
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert graph_counts == {
        "edges10_component6": 7,
        "edges11_component7": 15,
        "edges12_component8": 17,
    }
    assert maxima == {
        "edges10_component6": 250218343241809920,
        "edges11_component7": 253860219812966400,
        "edges12_component8": 255169391129395200,
    }
    assert all(value < q153 for value in maxima.values())
    return {
        "graph_type_counts": graph_counts,
        "allocation_maxima": allocation_maxima,
        "maxima": maxima,
        "maximum": max(maxima.values()),
        "over_q153": float(max(maxima.values()) / q153),
        "theorem": "All three tight active-order-six Q=126 triangle rows lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
