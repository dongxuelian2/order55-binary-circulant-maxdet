"""Exact support/Schur certificate placing Q=141 below Q=144."""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache

import networkx as nx

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import (
        SIZE13_Q90_MAXIMUM,
        weight_partitions,
    )
    from order15_mu3.scripts.verify_color_15_low_energy import (
        component_states,
        connected_types,
        q9_layer,
        weighted_layer,
    )
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q141_energy36_minimal_cross import (
        certificate as energy36_minimal_cross_certificate,
    )
    from order15_mu3.scripts.verify_size13_exact_loss import EXPECTED_LOCAL_LOSSES
    from order15_mu3.scripts.verify_support_independence import independence_number
    from order15_mu3.scripts.verify_trace_stability import sqrt_lower, stationary_upper
except ModuleNotFoundError:
    from verify_color_15_energy90 import SIZE13_Q90_MAXIMUM, weight_partitions
    from verify_color_15_low_energy import component_states, connected_types, q9_layer, weighted_layer
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q141_energy36_minimal_cross import certificate as energy36_minimal_cross_certificate
    from verify_size13_exact_loss import EXPECTED_LOCAL_LOSSES
    from verify_support_independence import independence_number
    from verify_trace_stability import sqrt_lower, stationary_upper


def _paw_energy36_upper(outside_norm: int) -> Fraction:
    """Phase-split residual bound for the unique four-edge paw support."""

    # For triangle gain z, put s=z+conj(z) in {2,-2,1,-1}.  The unit-weight
    # paw adjacency has characteristic polynomial
    #     p_s(t)=t^4-4t^2-s t+1.
    # The rational brackets below follow from p(a),p(b)>0, p'<0 to the left
    # of a, and p'>0 to the right of b.
    brackets = {
        2: (Fraction(-5, 3), Fraction(7, 3)),
        -2: (Fraction(-7, 3), Fraction(5, 3)),
        1: (Fraction(-11, 6), Fraction(7, 3)),
        -1: (Fraction(-13, 6), Fraction(2)),
    }

    def polynomial(t: Fraction, phase_sum: int) -> Fraction:
        return t**4 - 4 * t**2 - phase_sum * t + 1

    def derivative(t: Fraction, phase_sum: int) -> Fraction:
        return 4 * t**3 - 8 * t - phase_sum

    values = []
    for phase_sum, (lower, upper) in brackets.items():
        assert polynomial(lower, phase_sum) > 0 and derivative(lower, phase_sum) < 0
        assert polynomial(upper, phase_sum) > 0 and derivative(upper, phase_sum) > 0
        block_determinant = 81 * (526 + 5 * phase_sum) * 15**9
        lambda_min = 15 + 3 * lower
        lambda_max = 15 + 3 * upper
        trace_upper = Fraction(30) - Fraction(78) / lambda_max
        correction_upper = Fraction(39) / lambda_min
        residual_lower = sqrt_lower(Fraction(outside_norm)) - correction_upper
        assert residual_lower > 0
        values.append(block_determinant * ((trace_upper / 2) ** 2 - residual_lower**2))
    return max(values)


def _tree_type(graph: nx.Graph, norm9) -> dict[str, object]:
    """Return the component catalogue entry for an unlabeled q=9 tree."""

    graph = nx.convert_node_labels_to_integers(graph)
    edges = tuple(sorted(tuple(sorted(edge)) for edge in graph.edges()))
    states = component_states(graph, norm9)
    return {
        "vertices": graph.number_of_nodes(),
        "edges": graph.number_of_edges(),
        "edge_list": edges,
        "alpha": independence_number(graph.number_of_nodes(), edges),
        "has_leaf": True,
        "determinants": {state[0] for state in states},
        "states": states,
    }


def _energy18_upper() -> Fraction:
    loss = EXPECTED_LOCAL_LOSSES
    return max(
        3105 * 15**10 * (Fraction(15) - (Fraction(2) + loss["P3"][3105])) ** 2,
        216**2 * 15**9 * (
            Fraction(15) - (Fraction(9, 5) + 2 * loss["K2-q9"][216])
        ) ** 2,
    )


def _energy27_upper() -> Fraction:
    loss = EXPECTED_LOCAL_LOSSES
    values = [
        198 * 15**11 * (Fraction(15) - (Fraction(11, 5) + loss["K2-q27"][198])) ** 2,
        216**3 * 15**7 * (Fraction(15) - Fraction(12, 5)) ** 2,
        3105 * 216 * 15**8 * (
            Fraction(15) - (Fraction(8, 5) + loss["P3"][3105] + loss["K2-q9"][216])
        ) ** 2,
        44631 * 15**9 * (Fraction(15) - (Fraction(9, 5) + loss["P4"][44631])) ** 2,
        44550 * 15**9 * (Fraction(15) - (Fraction(9, 5) + loss["K1,3"][44550])) ** 2,
    ]
    values.extend(
        determinant_value * 15**10 * (Fraction(15) - (Fraction(2) + local_loss)) ** 2
        for determinant_value, local_loss in loss["K3"].items()
    )
    return max(values)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    q144 = max(stationary_upper(288, multiplicity) for multiplicity in range(1, 15))
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)

    # A size-13 class has internal energy 9,...,63.  Energies 9, 18, and 27
    # use the complete phase certificate and exact local inverse losses.
    # For larger energies, exact internal determinant layers (36,45) or the
    # support spectral caps (54,63) combine with the two-dimensional Schur
    # trace bound.  The five possible outside norms give exactly these rows.
    exact_internal = {}
    types5 = connected_types(5, norm9)
    for energy in (36, 45):
        exact_internal[energy] = max(
            weighted_layer(types5, weights, catalogue, order=13)["maximum_determinant"]
            for weights in weight_partitions(energy // 9)
        )
    types6 = connected_types(6, norm9)
    exact_internal[54] = max(
        weighted_layer(types6, weights, catalogue, order=13)["maximum_determinant"]
        for weights in weight_partitions(6)
    )
    types7 = connected_types(7, norm9)
    exact_internal[63] = max(
        weighted_layer(types7, weights, catalogue, order=13)["maximum_determinant"]
        for weights in weight_partitions(7)
    )
    pure36 = q9_layer(types5, 4, order=13)
    paw_profiles = [
        profile for profile in pure36["profiles"]
        if len(profile["component_edges"]) == 1
        and sorted(dict(nx.Graph(profile["component_edges"][0]).degree()).values()) == [1, 2, 2, 3]
    ]
    assert len(paw_profiles) == 1
    nonpaw36 = max(
        profile["maximum_determinant"]
        for profile in pure36["profiles"] if profile not in paw_profiles
    )
    assert nonpaw36 == 1674039150000000
    exact_tight36 = energy36_minimal_cross_certificate()
    remaining36 = max(
        profile["maximum_determinant"]
        for profile in pure36["profiles"]
        if profile not in paw_profiles and profile["maximum_determinant"] < nonpaw36
    )
    energy36_minimal_cross = max(
        Fraction(exact_tight36["maximum"]),
        Fraction(remaining36) * (Fraction(15) - Fraction(39, 21)) ** 2,
    )
    universal = {
        9: Fraction(SIZE13_Q90_MAXIMUM),
        18: _energy18_upper(),
        27: _energy27_upper(),
    }
    lambda_caps = {36: Fraction(22), 45: Fraction(23), 54: Fraction(24), 63: Fraction(25)}
    outside_norms = (0, 9, 27, 36, 63)
    size13_rows = []
    size13_values = []
    non_energy9_values = []
    for internal in range(9, 64, 9):
        for outside in outside_norms:
            cross = 141 - internal - outside
            if cross < 78 or cross % 9 != 6:
                continue
            if internal in universal:
                upper = universal[internal]
            else:
                block = (
                    Fraction(exact_internal[internal])
                    if internal in exact_internal
                    else stationary_upper_dimension(13, internal)
                )
                if internal == 36 and cross == 78 and outside == 27:
                    upper = energy36_minimal_cross
                else:
                    upper = block * (
                        Fraction(15) - Fraction(cross, 2) / lambda_caps[internal]
                    ) ** 2
            assert upper < q144
            size13_values.append(upper)
            if internal != 9:
                non_energy9_values.append(upper)
            size13_rows.append({
                "internal": internal,
                "cross": cross,
                "outside_norm": outside,
                "over_q144": float(upper / q144),
            })
    # This extra row is not needed at Q=141 itself, but makes the certificate
    # a reusable uniform bound for earlier color slices: Q=126 permits
    # internal energy 45 with the minimal cross energy 78.
    energy45_minimal_cross = Fraction(exact_internal[45]) * (
        Fraction(15) - Fraction(39, 23)
    ) ** 2
    assert energy45_minimal_cross < q144

    # If neither axis has a size-13 class, both color partitions are (14,1).
    # The equal-color-sign orbit lemma excludes an isolated internal support.
    # A no-isolate 99-energy support has at least seven edges, leaving exactly
    # four weight partitions.  The seven-edge case is a perfect matching.
    matching = 198**2 * 216**5
    layer8 = weighted_layer(
        connected_types(8, norm9), (4,) + (1,) * 7,
        catalogue, required_isolates=0, order=14,
    )
    layer9 = weighted_layer(
        connected_types(5, norm9), (3,) + (1,) * 8,
        catalogue, required_isolates=0, order=14,
    )
    # At eleven unit edges, connected_types(8) is complete except for a
    # possible nine-edge tree on ten vertices (followed by two K2's).
    types11 = list(connected_types(8, norm9))
    types11.extend(_tree_type(tree, norm9) for tree in nx.nonisomorphic_trees(10))
    layer11 = q9_layer(types11, 11, required_isolates=0, order=14)
    internal99 = {
        "3+3+1+1+1+1+1": matching,
        "4+1+1+1+1+1+1+1": layer8["maximum_determinant"],
        "3+1+1+1+1+1+1+1+1": layer9["maximum_determinant"],
        "eleven-unit": layer11["maximum_determinant"],
    }
    # Eleven support edges have clique number at most five.  Weighted
    # Motzkin--Straus gives rho(E)^2 <= 2*99*4/5 < (63/5)^2, so
    # lambda_max(G) < 138/5.
    lambda99 = Fraction(138, 5)
    assert Fraction(2 * 99 * 4, 5) < Fraction(63, 5) ** 2
    size14 = {
        name: Fraction(value) * (Fraction(15) - Fraction(42, 1) / lambda99)
        for name, value in internal99.items()
    }
    assert all(value < q144 for value in size14.values())
    certified_upper = max(*size13_values, energy45_minimal_cross, *size14.values())
    certified_excluding_energy9 = max(
        *non_energy9_values, energy45_minimal_cross, *size14.values()
    )

    return {
        "q144_envelope_numerator": q144.numerator,
        "q144_envelope_denominator": q144.denominator,
        "size13_rows": size13_rows,
        "size13_max_over_q144": max(row["over_q144"] for row in size13_rows),
        "size13_energy45_minimal_cross_over_q144": float(energy45_minimal_cross / q144),
        "size13_energy45_minimal_cross_numerator": energy45_minimal_cross.numerator,
        "size13_energy45_minimal_cross_denominator": energy45_minimal_cross.denominator,
        "size13_exact_internal_maxima": exact_internal,
        "size14_internal99_maxima": internal99,
        "size14_over_q144": {name: float(value / q144) for name, value in size14.items()},
        "certified_upper_numerator": certified_upper.numerator,
        "certified_upper_denominator": certified_upper.denominator,
        "certified_excluding_energy9_numerator": certified_excluding_energy9.numerator,
        "certified_excluding_energy9_denominator": certified_excluding_energy9.denominator,
        "theorem": "Every Q=141 Gram matrix lies below the Q=144 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
