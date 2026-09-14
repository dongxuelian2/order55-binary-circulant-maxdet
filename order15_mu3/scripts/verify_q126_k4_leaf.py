"""Exact phase certificate for the tight Q=126 weighted K4-plus-leaf branch."""

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
    catalogue = inner_product_values(15)
    graph = nx.complete_graph(4)
    graph.add_edge(0, 4)
    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    support_type = {
        "vertices": 5,
        "edges": 7,
        "edge_list": edges,
        "alpha": independence_number(5, edges),
        "has_leaf": True,
        "determinants": set(),
        "states": set(),
    }
    unit = weighted_layer(
        (support_type,), (1,) * 7, catalogue, required_isolates=0, order=5
    )["maximum_determinant"]
    weighted = weighted_layer(
        (support_type,), (3,) + (1,) * 6, catalogue, required_isolates=0, order=5
    )["maximum_determinant"]
    placements = {
        "heavy_outside": unit * 198 * 216**4,
        "heavy_in_component": weighted * 216**5,
    }
    maximum = max(placements.values())
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert unit == 594864 and weighted == 548451
    assert maximum == 257873424975691776 < q153
    return {
        "component_internal_maxima": {"unit": unit, "weighted": weighted},
        "full_placement_maxima": placements,
        "maximum": maximum,
        "over_q153": float(maximum / q153),
        "theorem": "The Q=126 K4-plus-leaf with five outside K2 components lies below Q=153.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
