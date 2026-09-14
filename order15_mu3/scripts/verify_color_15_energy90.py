"""Exact certificate for the all-same-color Q=90 boundary."""

from __future__ import annotations

import json
from functools import lru_cache

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_color_13_energy9 import arithmetic_admissible
    from order15_mu3.scripts.verify_color_15_energy81 import collect_kernel_states, rank_feasible_pairs
    from order15_mu3.scripts.verify_color_15_low_energy import connected_types, q9_layer, weighted_layer
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_bounds import is_rational_eisenstein_norm
    from verify_color_13_energy9 import arithmetic_admissible
    from verify_color_15_energy81 import collect_kernel_states, rank_feasible_pairs
    from verify_color_15_low_energy import connected_types, q9_layer, weighted_layer
    from verify_trace_stability import BENCHMARK


MIXED_PARTITIONS = (
    (9, 1),
    (7, 3),
    (7, 1, 1, 1),
    (4, 4, 1, 1),
    (4, 3, 3),
    (4, 3, 1, 1, 1),
    (4, 1, 1, 1, 1, 1, 1),
    (3, 3, 3, 1),
    (3, 3, 1, 1, 1, 1),
    (3, 1, 1, 1, 1, 1, 1, 1),
)
SIZE13_Q90_MAXIMUM = 289967495625000000


def weight_partitions(total: int, allowed=(9, 7, 4, 3, 1), maximum=None):
    """Descending integer partitions using the same-color norm units."""

    maximum = allowed[0] if maximum is None else maximum
    if total == 0:
        yield ()
        return
    for value in allowed:
        if value <= total and value <= maximum:
            for suffix in weight_partitions(total - value, allowed, value):
                yield (value,) + suffix


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    assert set(weight_partitions(10)) == {*MIXED_PARTITIONS, (1,) * 10}
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)

    # If either Gram support has an isolate, AH=HB supplies a mu3-valued
    # full-support kernel vector on every active component of the other
    # support.  Gauge that vector to one, enumerate the exact vertex-balance
    # equations, and then impose both trace/rank inequalities
    #     m_row v_col <= 15 d_col,  m_col v_row <= 15 d_row.
    isolated_layers = []
    mixed_counts = {}
    for weights in MIXED_PARTITIONS:
        layer = weighted_layer(
            connected_types(len(weights), norm9, kernel_only=True),
            weights,
            catalogue,
            mu3_kernel_only=True,
        )
        isolated_layers.append(layer)
        mixed_counts[str(weights)] = {
            "profiles": len(layer["profiles"]),
            "kernel_arithmetic_survivors": layer["intertwining_arithmetic_survivors"],
        }
    q9_isolated = q9_layer(connected_types(10, norm9, kernel_only=True), 10)
    isolated_layers.append(q9_isolated)
    states = collect_kernel_states(isolated_layers)
    rank_survivors = rank_feasible_pairs(states)
    assert states and not rank_survivors

    # With no isolates, ten unit edges are enumerated exactly.  If one edge
    # has norm 27 and seven have norm 9, degree sum 16 forces P3 + 6 K2.
    q9_no_isolate_layer = q9_layer(connected_types(5, norm9), 10)
    q9_no_isolate_profiles = [profile for profile in q9_no_isolate_layer["profiles"] if profile["isolates"] == 0]
    q9_no_isolate_values = sorted({
        value for profile in q9_no_isolate_profiles for value in profile["arithmetic_survivors"]
    })
    assert q9_no_isolate_values == [282438239494864896, 288608006862471168]
    mixed_no_isolate_values = sorted((2835 * 216**6, 3105 * 198 * 216**5))
    assert mixed_no_isolate_values == [287922477154959360, 289065026667479040]
    assert all(not is_rational_eisenstein_norm(value) for value in mixed_no_isolate_values)

    all_same_upper = max(q9_no_isolate_values)
    assert all_same_upper < SIZE13_Q90_MAXIMUM
    assert arithmetic_admissible(SIZE13_Q90_MAXIMUM)
    return {
        "energy": 90,
        "mixed_weight_partitions": MIXED_PARTITIONS,
        "isolated_mixed_counts": mixed_counts,
        "isolated_q9_profiles": len(q9_isolated["profiles"]),
        "isolated_kernel_states": {str(value): sorted(items) for value, items in sorted(states.items())},
        "isolated_rank_pair_survivors": rank_survivors,
        "no_isolate_q9_profiles": len(q9_no_isolate_profiles),
        "no_isolate_q9_arithmetic_values": q9_no_isolate_values,
        "no_isolate_mixed_values": mixed_no_isolate_values,
        "all_same_color_upper": all_same_upper,
        "color_13_1_1_upper": SIZE13_Q90_MAXIMUM,
        "all_Q90_upper": SIZE13_Q90_MAXIMUM,
        "benchmark": BENCHMARK,
        "theorem": "Every order-15 Gram matrix at total energy Q=90 has determinant at most 289967495625000000.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
