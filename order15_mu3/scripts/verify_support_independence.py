"""Replay support-graph consequences of B_3(15,10)=12."""

from __future__ import annotations

import json
from itertools import combinations


CODE_BOUND = 12
SOURCE_DOI = "10.28919/jmcs/4964"


def independence_number(vertices: int, edges: tuple[tuple[int, int], ...]) -> int:
    edge_sets = {frozenset(edge) for edge in edges}
    for size in range(vertices, -1, -1):
        for subset in combinations(range(vertices), size):
            if all(frozenset(pair) not in edge_sets for pair in combinations(subset, 2)):
                return size
    raise AssertionError("unreachable")


def main() -> None:
    profiles = {
        "energy9_one_edge": (14, ((0, 1),)),
        "energy18_disjoint": (14, ((0, 1), (2, 3))),
        "energy18_path": (14, ((0, 1), (1, 2))),
        "energy27_one_edge": (14, ((0, 1),)),
        "energy27_matching3": (14, ((0, 1), (2, 3), (4, 5))),
        "energy27_path_plus_edge": (14, ((0, 1), (1, 2), (3, 4))),
        "energy27_path4": (14, ((0, 1), (1, 2), (2, 3))),
        "energy27_star4": (14, ((0, 1), (0, 2), (0, 3))),
        "energy27_triangle": (14, ((0, 1), (1, 2), (0, 2))),
    }
    alphas = {name: independence_number(*profile) for name, profile in profiles.items()}
    assert alphas == {
        "energy9_one_edge": 13,
        "energy18_disjoint": 12,
        "energy18_path": 13,
        "energy27_one_edge": 13,
        "energy27_matching3": 11,
        "energy27_path_plus_edge": 12,
        "energy27_path4": 12,
        "energy27_star4": 13,
        "energy27_triangle": 12,
    }
    print(json.dumps({
        "external_exact_code_bound": "B_3(15,10)=12",
        "source_doi": SOURCE_DOI,
        "support_independence_numbers": alphas,
        "excluded_14_class_profiles": [name for name, alpha in alphas.items() if alpha > CODE_BOUND],
        "minimum_internal_energies_by_class_size": {"13": 9, "14": 18, "15": 27},
    }, indent=2))


if __name__ == "__main__":
    main()
