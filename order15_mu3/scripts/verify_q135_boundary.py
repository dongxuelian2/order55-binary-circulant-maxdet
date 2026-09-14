"""Exact spectral/support certificate placing Q=135 below Q=141."""

from __future__ import annotations

import json
from fractions import Fraction

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import SIZE13_Q90_MAXIMUM, weight_partitions
    from order15_mu3.scripts.verify_color_15_low_energy import connected_types, weighted_layer
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q126_boundary import capped_trace_candidates
    from order15_mu3.scripts.verify_size13_exact_loss import EXPECTED_LOCAL_LOSSES
    from order15_mu3.scripts.verify_trace_stability import sqrt_lower, stationary_upper
except ModuleNotFoundError:
    from verify_color_15_energy90 import SIZE13_Q90_MAXIMUM, weight_partitions
    from verify_color_15_low_energy import connected_types, weighted_layer
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q126_boundary import capped_trace_candidates
    from verify_size13_exact_loss import EXPECTED_LOCAL_LOSSES
    from verify_trace_stability import sqrt_lower, stationary_upper


def fixed_isolate_capped_upper(isolates: int, variance_energy: int, cap: Fraction) -> Fraction:
    """KKT product upper with `isolates` eigenvalues fixed at 15."""

    order = 15 - isolates
    trace = Fraction(15 * order)
    trace_square = Fraction(225 * order + variance_energy)
    candidates = []
    for multiplicity in range(1, order):
        t = sqrt_lower(Fraction(variance_energy, order * multiplicity * (order - multiplicity)))
        high = Fraction(15) + (order - multiplicity) * t
        low = Fraction(15) - multiplicity * t
        if high <= cap and low > 0:
            candidates.append(high**multiplicity * low ** (order - multiplicity))
    for capped in range(1, order):
        remaining = order - capped
        remaining_trace = trace - capped * cap
        remaining_square = trace_square - capped * cap**2
        variance = remaining_square - remaining_trace**2 / remaining
        if variance < 0:
            continue
        for multiplicity in range(1, remaining):
            t = sqrt_lower(Fraction(variance, remaining * multiplicity * (remaining - multiplicity)))
            high = remaining_trace / remaining + (remaining - multiplicity) * t
            low = remaining_trace / remaining - multiplicity * t
            if high <= cap and low > 0:
                candidates.append(cap**capped * high**multiplicity * low ** (remaining - multiplicity))
    assert candidates
    return 15**isolates * max(candidates)


def certificate() -> dict[str, object]:
    q141 = max(stationary_upper(282, multiplicity) for multiplicity in range(1, 15))

    # All-same Q=135 has at most 15 support edges and clique number at most
    # six.  The weighted Motzkin--Straus bound is rho(E)^2<=225, hence
    # lambda_max(G)<=30.  Exact KKT boundary enumeration closes the case.
    all_same = max(row[3] for row in capped_trace_candidates(Fraction(30), 270))
    assert all_same < q141

    # For (13,1,1), exact energy-9 phases and exact sparse internal
    # determinants suffice once the cross energy is bounded below by 78.
    assert SIZE13_Q90_MAXIMUM < q141
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    internal_types = connected_types(5, norm9)
    size13 = {"energy9": Fraction(SIZE13_Q90_MAXIMUM)}
    loss = EXPECTED_LOCAL_LOSSES
    size13["energy18"] = max(
        3105 * 15**10 * (Fraction(15) - Fraction(57, 23)) ** 2,
        216**2 * 15**9 * (Fraction(15) - Fraction(37, 15)) ** 2,
    )
    energy27_exact = [
        198 * 15**11 * (Fraction(15) - (Fraction(11, 5) + loss["K2-q27"][198])) ** 2,
        216**3 * 15**7 * (Fraction(15) - Fraction(12, 5)) ** 2,
        3105 * 216 * 15**8 * (
            Fraction(15) - (Fraction(8, 5) + loss["P3"][3105] + loss["K2-q9"][216])
        ) ** 2,
        44631 * 15**9 * (Fraction(15) - (Fraction(9, 5) + loss["P4"][44631])) ** 2,
        44550 * 15**9 * (Fraction(15) - (Fraction(9, 5) + loss["K1,3"][44550])) ** 2,
    ]
    energy27_exact += [
        determinant_value * 15**10 * (
            Fraction(15) - (Fraction(2) + local_loss)
        ) ** 2
        for determinant_value, local_loss in loss["K3"].items()
    ]
    size13["energy27"] = max(energy27_exact)
    lambda_by_energy = {36: Fraction(22), 45: Fraction(23)}
    for energy, lambda_upper in lambda_by_energy.items():
        maximum = max(
            weighted_layer(internal_types, weights, catalogue, order=13)["maximum_determinant"]
            for weights in weight_partitions(energy // 9)
        )
        size13[f"energy{energy}"] = maximum * (
            Fraction(15) - Fraction(39, 1) / lambda_upper
        ) ** 2
    size13["energy54"] = stationary_upper_dimension(13, 54) * (
        Fraction(15) - Fraction(39, 24)
    ) ** 2
    assert all(value < q141 for value in size13.values())

    return {
        "q141_envelope_numerator": q141.numerator,
        "q141_envelope_denominator": q141.denominator,
        "all_same_lambda_cap": 30,
        "all_same_over_q141": float(all_same / q141),
        "size13_over_q141": {key: float(value / q141) for key, value in size13.items()},
        "theorem": "Every Q=135 Gram matrix lies below the Q=141 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
