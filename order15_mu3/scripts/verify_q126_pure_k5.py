"""Exact phase certificate for the pure Q=126 isolated K5 component."""

from __future__ import annotations

import json
from functools import lru_cache

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


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    graph = nx.complete_graph(5)
    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    support_type = {
        "vertices": 5,
        "edges": 10,
        "edge_list": edges,
        "alpha": independence_number(5, edges),
        "has_leaf": False,
        "determinants": set(),
        "states": set(),
    }
    layer = weighted_layer(
        (support_type,),
        (1,) * 10,
        inner_product_values(15),
        required_isolates=0,
        order=5,
    )
    component_maximum = layer["maximum_determinant"]
    joined0 = component_maximum * 216**4 * 15**2
    leaf_graph = nx.complete_graph(5)
    leaf_graph.add_edge(0, 5)
    leaf_edges = tuple(sorted(tuple(sorted(edge)) for edge in leaf_graph.edges()))
    leaf_type = {
        "vertices": 6,
        "edges": 11,
        "edge_list": leaf_edges,
        "alpha": independence_number(6, leaf_edges),
        "has_leaf": True,
        "determinants": set(),
        "states": set(),
    }
    leaf_layer = weighted_layer(
        (leaf_type,),
        (1,) * 11,
        inner_product_values(15),
        required_isolates=0,
        order=6,
    )
    leaf_component_maximum = leaf_layer["maximum_determinant"]
    maxima = {
        "joined0": joined0,
        "joined1_leaf": leaf_component_maximum * 216**3 * 15**3,
    }
    maximum = max(maxima.values())
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert component_maximum == 559872
    assert leaf_component_maximum == 8024832
    assert joined0 == 274211883004723200
    assert maxima["joined1_leaf"] == 272942383546368000
    assert maximum < q153
    return {
        "phase_profile_count": len(layer["profiles"]),
        "component_maximum": component_maximum,
        "leaf_component_maximum": leaf_component_maximum,
        "maxima": maxima,
        "maximum": maximum,
        "over_q153": float(maximum / q153),
        "theorem": "The pure Q=126 isolated K5 and K5-plus-leaf branches lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
