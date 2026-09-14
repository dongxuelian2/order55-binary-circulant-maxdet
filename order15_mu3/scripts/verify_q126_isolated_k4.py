"""Exact phase certificate for the isolated-K4 Q=126 branches."""

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
    graph = nx.complete_graph(4)
    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    support_type = {
        "vertices": 4,
        "edges": 6,
        "edge_list": edges,
        "alpha": independence_number(4, edges),
        "has_leaf": False,
        "determinants": set(),
        "states": set(),
    }
    catalogue = inner_product_values(15)
    layers = {
        "unit": weighted_layer(
            (support_type,), (1,) * 6, catalogue, required_isolates=0, order=4
        ),
        "weighted3": weighted_layer(
            (support_type,), (3,) + (1,) * 5, catalogue, required_isolates=0, order=4
        ),
        "weighted4": weighted_layer(
            (support_type,), (4,) + (1,) * 5, catalogue, required_isolates=0, order=4
        ),
    }
    path = nx.path_graph(3)
    path_edges = tuple(sorted(tuple(sorted(edge)) for edge in path.edges()))
    path_type = {
        "vertices": 3,
        "edges": 2,
        "edge_list": path_edges,
        "alpha": independence_number(3, path_edges),
        "has_leaf": True,
        "determinants": set(),
        "states": set(),
    }
    path_unit = weighted_layer(
        (path_type,), (1, 1), catalogue, required_isolates=0, order=3
    )["maximum_determinant"]
    path_weighted3 = weighted_layer(
        (path_type,), (3, 1), catalogue, required_isolates=0, order=3
    )["maximum_determinant"]
    component_maxima = {
        name: layer["maximum_determinant"] for name, layer in layers.items()
    }
    maxima = {
        "edges11_weight4": max(
            component_maxima["unit"] * 189 * 216**4 * 15,
            component_maxima["weighted4"] * 216**5 * 15,
        ),
        "edges12_weight3": max(
            component_maxima["weighted3"] * path_unit * 216**4,
            component_maxima["unit"] * path_weighted3 * 216**4,
            component_maxima["unit"] * path_unit * 198 * 216**3,
        ),
    }
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert component_maxima["unit"] == 41472
    assert component_maxima["weighted3"] == 38232
    assert component_maxima["weighted4"] == 36936
    assert maxima["edges11_weight4"] == 260501288854487040
    assert maxima["edges12_weight3"] == 258406614748200960
    assert all(value < q153 for value in maxima.values())
    return {
        "component_maxima": component_maxima,
        "maxima": maxima,
        "maximum": max(maxima.values()),
        "over_q153": float(max(maxima.values()) / q153),
        "theorem": "The isolated-K4 Q=126 branches, including the forced P3 case, lie below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
