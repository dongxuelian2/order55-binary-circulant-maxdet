"""Exact determinant audit for the tight Q=126 triangle components."""

from __future__ import annotations

import json
from functools import lru_cache

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_low_energy import (
        connected_near_tree_graphs,
        connected_unicyclic_graphs,
        q9_layer,
        weighted_layer,
    )
    from order15_mu3.scripts.verify_support_independence import independence_number
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_color_15_low_energy import connected_near_tree_graphs, connected_unicyclic_graphs, q9_layer, weighted_layer
    from verify_support_independence import independence_number
    from verify_trace_stability import stationary_upper


def _type(graph: nx.Graph) -> dict[str, object]:
    graph = nx.convert_node_labels_to_integers(graph)
    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    return {
        "vertices": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "edge_list": edges,
        "alpha": independence_number(graph.number_of_nodes(), edges),
        "has_leaf": any(graph.degree(vertex) == 1 for vertex in graph),
        "determinants": set(),
        "states": set(),
    }


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    catalogue = inner_product_values(15)
    unicyclic = tuple(
        graph for graph in connected_unicyclic_graphs(9)
        if len(nx.cycle_basis(graph)[0]) == 3
    )
    unicyclic10 = tuple(
        graph for graph in connected_unicyclic_graphs(10)
        if len(nx.cycle_basis(graph)[0]) == 3
    )
    bicyclic = tuple(
        graph for graph in connected_near_tree_graphs(8, 9)
        if any(len(cycle) == 3 for cycle in nx.cycle_basis(graph))
        and max(map(len, nx.find_cliques(graph))) < 4
    )
    bicyclic7 = tuple(
        graph for graph in connected_near_tree_graphs(7, 8)
        if any(len(cycle) == 3 for cycle in nx.cycle_basis(graph))
        and max(map(len, nx.find_cliques(graph))) < 4
    )

    def layers_for(graphs, order, edge_count):
        types = tuple(_type(graph) for graph in graphs)
        return {
            "unit_component": weighted_layer(
                types,
                (1,) * edge_count,
                catalogue,
                required_isolates=0,
                order=order,
            ),
            "weighted_component": weighted_layer(
                types,
                (3,) + (1,) * (edge_count - 1),
                catalogue,
                required_isolates=0,
                order=order,
            ),
        }

    layers = {
        "active10_unicyclic": layers_for(unicyclic10, 10, 10),
        "active9_unicyclic": layers_for(unicyclic, 9, 9),
        "active8_bicyclic": layers_for(bicyclic, 8, 9),
        "active7_bicyclic": layers_for(bicyclic7, 7, 8),
    }
    component_edges = {
        "active10_unicyclic": 10,
        "active9_unicyclic": 9,
        "active8_bicyclic": 9,
        "active7_bicyclic": 8,
    }
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    maxima = {}
    for name, pair in layers.items():
        edges = component_edges[name]
        outside_edges = 12 - edges
        isolates = 15 - int(name.removeprefix("active").split("_")[0]) - 2 * outside_edges
        assert isolates >= 0
        maxima[name] = max(
            pair["unit_component"]["maximum_determinant"]
            * 198
            * 216 ** (outside_edges - 1)
            * 15**isolates,
            pair["weighted_component"]["maximum_determinant"]
            * 216**outside_edges
            * 15**isolates,
        )
    assert len(unicyclic10) == 299
    assert len(unicyclic) == 117 and len(bicyclic) == 173 and len(bicyclic7) == 51
    assert maxima == {
        "active10_unicyclic": 246462947301769920,
        "active9_unicyclic": 246911098793504256,
        "active8_bicyclic": 251108976455562240,
        "active7_bicyclic": 251565599541989376,
    }
    assert all(value < q153 for value in maxima.values())
    return {
        "graph_type_counts": {
            "active10_unicyclic": len(unicyclic10),
            "active9_unicyclic": len(unicyclic),
            "active8_bicyclic": len(bicyclic),
            "active7_bicyclic": len(bicyclic7),
        },
        "profile_counts": {
            name: {kind: len(layer["profiles"]) for kind, layer in pair.items()}
            for name, pair in layers.items()
        },
        "maxima": maxima,
        "maximum": max(maxima.values()),
        "over_q153": float(max(maxima.values()) / q153),
        "theorem": "The tight Q=126 K4-free triangle-component branches lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
