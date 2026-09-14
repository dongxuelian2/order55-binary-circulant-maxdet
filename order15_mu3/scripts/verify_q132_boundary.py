"""Exact determinant-matching certificate placing Q=132 below Q=135."""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from order15_mu3.scripts.verify_color_15_low_energy import connected_types, q9_layer, weighted_layer
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_boundary import generic_size13
    from order15_mu3.scripts.verify_q114_internal_support import EXPECTED_MAXIMA
    from order15_mu3.scripts.verify_q117_boundary import schur_two_by_two
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, sqrt_lower, stationary_upper
except ModuleNotFoundError:
    from verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from verify_color_15_low_energy import connected_types, q9_layer, weighted_layer
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_boundary import generic_size13
    from verify_q114_internal_support import EXPECTED_MAXIMA
    from verify_q117_boundary import schur_two_by_two
    from verify_trace_stability import BENCHMARK, sqrt_lower, stationary_upper


def schur_with_residual(
    block_determinant: int,
    lambda_max_upper: Fraction,
    lambda_min_lower: Fraction,
    outside_norm: int,
) -> Fraction:
    """Two-hub Schur bound when all 26 cross norms are minimal."""

    trace_upper = Fraction(30) - Fraction(78, 1) / lambda_max_upper
    correction_upper = Fraction(39, 1) / lambda_min_lower
    residual_lower = sqrt_lower(Fraction(outside_norm)) - correction_upper
    assert residual_lower > 0
    return block_determinant * ((trace_upper / 2) ** 2 - residual_lower**2)


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    q135 = max(stationary_upper(270, multiplicity) for multiplicity in range(1, 15))

    # Exact size-13 refinements.  Energy 9 is covered by the prior complete
    # phase enumeration.  At (e,c,r)=(18,78,36), retain the outside Schur
    # residual.  At (27,78,27), all five forest shapes use lambda_min>=9;
    # triangle phases are split by their exact determinant/spectral classes.
    assert SIZE13_Q90_MAXIMUM < q135
    lambda_p3 = Fraction(193, 10)
    energy18 = {
        "cross87": max(
            schur_two_by_two(3105 * 15**10, lambda_p3, 87),
            schur_two_by_two(216**2 * 15**9, Fraction(18), 87),
        ),
        "cross78-outside36-P3": schur_with_residual(
            3105 * 15**10, lambda_p3, Fraction(107, 10), 36
        ),
        "cross78-outside36-2K2": schur_with_residual(
            216**2 * 15**9, Fraction(18), Fraction(12), 36
        ),
    }
    forest27 = {
        "norm27-K2": (198 * 15**11, Fraction(21)),
        "3K2": (216**3 * 15**7, Fraction(18)),
        "P3+K2": (3105 * 216 * 15**8, lambda_p3),
        "P4": (44631 * 15**9, Fraction(20)),
        "K1,3": (44550 * 15**9, Fraction(21)),
    }
    energy27 = {
        name: schur_with_residual(det, lmax, Fraction(9), 27)
        for name, (det, lmax) in forest27.items()
    }
    triangle_classes = (
        (2916, Fraction(9)), (2943, Fraction(9)),
        (2997, Fraction(10)), (3024, Fraction(12)),
    )
    # For determinant 2997, p(10)=-17 and p'>0 on [9,10], proving
    # lambda_min>10.  The 3024 polynomial factors as (x-12)^2(x-21).
    assert 10**3 - 45 * 10**2 + 648 * 10 - 2997 == -17
    for determinant_value, lmin in triangle_classes:
        energy27[f"K3-det{determinant_value}"] = schur_with_residual(
            determinant_value * 15**10, Fraction(21), lmin, 27
        )
    assert all(value < q135 for value in energy18.values())
    assert all(value < q135 for value in energy27.values())

    # Every other size-13 allocation is handled by the generic bounds.
    outside_norms = (0, 9, 27, 36)
    generic_rows = []
    for internal in range(9, 55, 9):
        for outside in outside_norms:
            cross = 132 - internal - outside
            if cross < 78 or cross % 9 != 6:
                continue
            if internal == 9 or (internal == 18 and cross <= 87) or (internal == 27 and cross == 78):
                continue
            if internal <= 36:
                upper = generic_size13(internal, cross)
            else:
                lambda_upper = Fraction(23 if internal == 45 else 24)
                upper = stationary_upper_dimension(13, internal) * (
                    Fraction(15) - Fraction(cross, 2) / lambda_upper
                ) ** 2
            assert upper < q135
            generic_rows.append({
                "internal": internal, "cross": cross, "outside_norm": outside,
                "over_q135": float(upper / q135),
            })

    # A determinant above q135 cannot have a size-13 partition on either
    # axis.  Thus both axes are (14,1).  Internally isolated supports are
    # impossible by the equal-color-sign orbit lemma.  The no-isolate cases
    # at e=63,72,81,90 reduce to friendship supports or the exact layers here.
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    e81 = q9_layer(connected_types(5, norm9), 9, required_isolates=0, order=14)
    e90_len8 = weighted_layer(
        connected_types(8, norm9), (3,) + (1,) * 7,
        catalogue, required_isolates=0, order=14,
    )
    e90_len10 = q9_layer(connected_types(7, norm9), 10, required_isolates=0, order=14)
    assert e81["maximum_determinant"] == 20438941279518720
    assert e90_len8["maximum_determinant"] == 19237545177523200
    assert e90_len10["maximum_determinant"] == 19905751507009536
    no_isolate = {
        "friendship": Fraction(BENCHMARK),
        "e72": EXPECTED_MAXIMA[72] * (Fraction(15) - Fraction(60, 1) / Fraction(127, 5)),
        "e81-unit9": e81["maximum_determinant"] * (Fraction(15) - Fraction(51, 1) / Fraction(261, 10)),
        "e90-weighted8": e90_len8["maximum_determinant"] * (Fraction(15) - Fraction(42, 27)),
        "e90-unit10": e90_len10["maximum_determinant"] * (Fraction(15) - Fraction(42, 27)),
    }
    assert all(value < q135 for value in no_isolate.values())

    return {
        "q135_envelope_numerator": q135.numerator,
        "q135_envelope_denominator": q135.denominator,
        "size13_energy9_maximum": SIZE13_Q90_MAXIMUM,
        "size13_energy18_over_q135": {key: float(value / q135) for key, value in energy18.items()},
        "size13_energy27_over_q135": {key: float(value / q135) for key, value in energy27.items()},
        "size13_generic": generic_rows,
        "size14_no_isolate_over_q135": {key: float(value / q135) for key, value in no_isolate.items()},
        "theorem": "Every Q=132 Gram matrix lies below the Q=135 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
