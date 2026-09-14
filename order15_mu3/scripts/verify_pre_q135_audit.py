"""Compact exact audit of every nonempty energy below Q=135."""

from __future__ import annotations

import json
from fractions import Fraction

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from order15_mu3.scripts.verify_color_15_low_energy import connected_types, weighted_layer
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_internal_support import EXPECTED_MAXIMA
    from order15_mu3.scripts.verify_q126_boundary import capped_trace_candidates
    from order15_mu3.scripts.verify_size13_exact_loss import certificate as exact_loss_certificate
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, stationary_upper
except ModuleNotFoundError:
    from verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from verify_color_15_low_energy import connected_types, weighted_layer
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_internal_support import EXPECTED_MAXIMA
    from verify_q126_boundary import capped_trace_candidates
    from verify_size13_exact_loss import certificate as exact_loss_certificate
    from verify_trace_stability import BENCHMARK, stationary_upper


def certificate() -> dict[str, object]:
    q135 = max(stationary_upper(270, multiplicity) for multiplicity in range(1, 15))

    # Size-13 partitions.  Energy 9 has a complete phase-enumeration cap;
    # minimal-cross energies 18 and 27 use exact local quadratic losses.
    losses = exact_loss_certificate()
    assert SIZE13_Q90_MAXIMUM < q135
    assert max(losses["energy18_bounds_over_q135"].values()) < 1
    assert max(losses["energy27_bounds_over_q135"].values()) < 1

    # At energy 36, exact internal determinants replace the loose stationary
    # determinant.  The weighted spectral radii are at most 6,6,<7 for the
    # three partitions, giving lambda maxima 21,21,22.
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    types4 = connected_types(4, norm9)
    energy36 = {}
    for weights, lambda_upper in (
        ((4,), Fraction(21)), ((3, 1), Fraction(21)), ((1, 1, 1, 1), Fraction(22)),
    ):
        internal_maximum = weighted_layer(types4, weights, catalogue, order=13)["maximum_determinant"]
        upper = internal_maximum * (Fraction(15) - Fraction(39, 1) / lambda_upper) ** 2
        assert upper < q135
        energy36[str(weights)] = {
            "internal_maximum": internal_maximum,
            "over_q135": float(upper / q135),
        }
    # Higher internal energies arise only with cross energy at least 78; the
    # weighted clique bounds lambda<23 and lambda<24 suffice.
    high_internal = {
        45: stationary_upper_dimension(13, 45) * (Fraction(15) - Fraction(39, 23)) ** 2,
        54: stationary_upper_dimension(13, 54) * (Fraction(15) - Fraction(39, 24)) ** 2,
    }
    assert all(value < q135 for value in high_internal.values())

    # All-same Q=126: reuse the exact KKT list under lambda_max<=146/5.
    q126_all_same = max(row[3] for row in capped_trace_candidates(Fraction(146, 5)))
    assert q126_all_same < q135

    # Size-14 determinant matching: once every size-13 alternative is below
    # q135, a larger value forces (14,1) on both axes.  Internally isolated
    # pairs contradict the equal-color-sign orbit lemma.  The remaining
    # no-isolate supports before Q=132 are friendship, energy-72, or the
    # nine-unit-edge energy-81 layer.
    e72 = EXPECTED_MAXIMA[72] * (Fraction(15) - Fraction(42, 1) / Fraction(127, 5))
    e81_layer = __import__(
        "order15_mu3.scripts.verify_color_15_low_energy", fromlist=["q9_layer"]
    ).q9_layer(connected_types(5, norm9), 9, required_isolates=0, order=14)
    e81 = e81_layer["maximum_determinant"] * (
        Fraction(15) - Fraction(42, 1) / Fraction(261, 10)
    )
    assert e81_layer["maximum_determinant"] == 20438941279518720
    assert BENCHMARK < q135 and e72 < q135 and e81 < q135

    return {
        "q135_envelope_numerator": q135.numerator,
        "q135_envelope_denominator": q135.denominator,
        "size13_energy9_maximum": SIZE13_Q90_MAXIMUM,
        "size13_energy36": energy36,
        "size13_high_internal_over_q135": {str(key): float(value / q135) for key, value in high_internal.items()},
        "q126_all_same_over_q135": float(q126_all_same / q135),
        "size14_no_isolate_over_q135": {
            "friendship": float(Fraction(BENCHMARK) / q135),
            "energy72": float(e72 / q135),
            "energy81": float(e81 / q135),
        },
        "theorem": "Every strict counterexample with Q<135 lies below the Q=135 envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
