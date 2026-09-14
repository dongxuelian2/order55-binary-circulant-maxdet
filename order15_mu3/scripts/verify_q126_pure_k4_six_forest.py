"""Exact phase certificate for pure Q=126 K4 plus a six-edge rooted forest."""

from __future__ import annotations

import json
from functools import lru_cache
from itertools import product

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
    """Enumerate by contracting the distinguished K4 to a rooted tree vertex."""

    buckets: dict[str, list[nx.Graph]] = {}
    for tree in nx.nonisomorphic_trees(7):
        for root in tree.nodes():
            neighbours = tuple(tree.neighbors(root))
            outside = tuple(vertex for vertex in tree.nodes() if vertex != root)
            relabel = {vertex: index + 4 for index, vertex in enumerate(outside)}
            for attachments in product(range(4), repeat=len(neighbours)):
                graph = nx.complete_graph(4)
                graph.add_nodes_from(range(4, 10))
                attachment = dict(zip(neighbours, attachments))
                for left, right in tree.edges():
                    if left == root:
                        graph.add_edge(attachment[right], relabel[right])
                    elif right == root:
                        graph.add_edge(attachment[left], relabel[left])
                    else:
                        graph.add_edge(relabel[left], relabel[right])
                key = nx.weisfeiler_lehman_graph_hash(graph)
                bucket = buckets.setdefault(key, [])
                if not any(nx.is_isomorphic(graph, prior) for prior in bucket):
                    bucket.append(graph)
    result = []
    for graph in (item for bucket in buckets.values() for item in bucket):
        edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
        assert graph.number_of_nodes() == 10 and graph.number_of_edges() == 12
        result.append({
            "vertices": 10,
            "edges": 12,
            "edge_list": edges,
            "alpha": independence_number(10, edges),
            "has_leaf": any(graph.degree(vertex) == 1 for vertex in graph),
            "determinants": set(),
            "states": set(),
        })
    return tuple(result)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    types = _support_types()
    layer = weighted_layer(
        types,
        (1,) * 12,
        inner_product_values(15),
        required_isolates=0,
        order=10,
    )
    component_maximum = layer["maximum_determinant"]
    maximum = component_maximum * 216**2 * 15
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert len(types) == 124
    assert component_maximum == 364958721792
    assert maximum == 255412711858913280
    assert maximum < q153
    return {
        "graph_type_count": len(types),
        "profile_count": len(layer["profiles"]),
        "component_maximum": component_maximum,
        "maximum": maximum,
        "over_q153": float(maximum / q153),
        "theorem": "The pure Q=126 K4 plus six-edge rooted-forest branch lies below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
