"""Exact Schur certificate placing every Q=123 color case below Q=126."""

from __future__ import annotations

import json
from fractions import Fraction

from maxdet.mu3 import inner_product_values

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from order15_mu3.scripts.verify_color_15_low_energy import connected_types, q9_layer
    from order15_mu3.scripts.verify_q114_internal_support import EXPECTED_MAXIMA
    from order15_mu3.scripts.verify_q114_boundary import LAMBDA_14, generic_size13, one_edge_tight
    from order15_mu3.scripts.verify_q117_boundary import schur_two_by_two
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from verify_color_15_low_energy import connected_types, q9_layer
    from verify_q114_internal_support import EXPECTED_MAXIMA
    from verify_q114_boundary import LAMBDA_14, generic_size13, one_edge_tight
    from verify_q117_boundary import schur_two_by_two
    from verify_trace_stability import stationary_upper


def certificate() -> dict[str, object]:
    q126 = max(stationary_upper(252, multiplicity) for multiplicity in range(1, 15))
    q132 = max(stationary_upper(264, multiplicity) for multiplicity in range(1, 15))

    rows14 = []
    lambda14 = {9: Fraction(18), **LAMBDA_14, 81: Fraction(261, 10)}
    for internal in range(9, 82, 9):
        cross = 123 - internal
        bound = lambda14[internal]
        if internal == 81:
            # There are at most nine nonzero internal edges, hence no K5 and
            # clique number <=4.  Weighted Cauchy--Schwarz plus
            # Motzkin--Straus gives rho(E)^2<=2*81*(1-1/4)=243/2.
            assert Fraction(243, 2) < (bound - 15) ** 2
        upper = stationary_upper_dimension(14, internal) * (Fraction(15) - Fraction(cross, 1) / bound)
        assert upper < q126
        if internal <= 63:
            assert upper < q132
        rows14.append({
            "internal": internal,
            "cross": cross,
            "lambda_upper": [bound.numerator, bound.denominator],
            "over_q126": float(upper / q126),
        })

    # For (13,2), the two-row block's attainable norms through 36 are
    # 0,9,27,36.  These and c>=78 leave exactly the following allocations.
    allocations13 = []
    for internal in range(9, 46, 9):
        for outside in (0, 9, 27, 36):
            cross = 123 - internal - outside
            if cross >= 78:
                allocations13.append((internal, cross, outside))
    assert len(allocations13) == 12

    tight = {(9, 78, 36), (18, 78, 27), (36, 78, 9)}
    generic_rows = []
    for internal, cross, outside in allocations13:
        if (internal, cross, outside) in tight:
            continue
        if internal <= 36:
            upper = generic_size13(internal, cross)
        else:
            # e=45 has at most five edges and hence clique number at most
            # three.  Thus rho(E)^2<=2*45*(1-1/3)=60<8^2.
            assert (Fraction(23) - 15) ** 2 > 60
            upper = stationary_upper_dimension(13, 45) * (Fraction(15) - Fraction(39, 23)) ** 2
        assert upper < q126
        generic_rows.append({
            "internal": internal,
            "cross": cross,
            "outside_norm": outside,
            "over_q126": float(upper / q126),
        })

    tight_bounds = {
        (9, 78, 36): one_edge_tight(36),
        # At e=18 the only supports are P3 and 2K2.
        (18, 78, 27): max(
            schur_two_by_two(3105 * 15**10, Fraction(193, 10)),
            schur_two_by_two(216**2 * 15**9, Fraction(18)),
        ),
        # At e=36, at most four internal edges and clique number <=3 imply
        # rho(E)^2<=48<7^2 and lambda_max(A)<22.
        (36, 78, 9): stationary_upper_dimension(13, 36) * (Fraction(15) - Fraction(39, 22)) ** 2,
    }
    assert all(value < q126 for value in tight_bounds.values())
    tight_rows = [
        {
            "internal": allocation[0],
            "cross": allocation[1],
            "outside_norm": allocation[2],
            "over_q126": float(value / q126),
        }
        for allocation, value in tight_bounds.items()
    ]

    # Stronger audit for advancing the global envelope to Q=132.  Every
    # size-13 energy-9 case is bounded by its exact phase enumeration.  The
    # remaining tight internal supports use the same exact shape bounds.
    assert SIZE13_Q90_MAXIMUM < q132
    lambda_p3 = Fraction(193, 10)
    size13_q132 = {
        "energy9": Fraction(SIZE13_Q90_MAXIMUM),
        "energy18": max(
            schur_two_by_two(3105 * 15**10, lambda_p3),
            schur_two_by_two(216**2 * 15**9, Fraction(18)),
        ),
        "energy27": max(
            schur_two_by_two(198 * 15**11, Fraction(21)),
            schur_two_by_two(216**3 * 15**7, Fraction(18)),
            schur_two_by_two(3105 * 216 * 15**8, lambda_p3),
            schur_two_by_two(44631 * 15**9, Fraction(20)),
            schur_two_by_two(44550 * 15**9, Fraction(21)),
            schur_two_by_two(3024 * 15**10, Fraction(21)),
        ),
        "energy36": stationary_upper_dimension(13, 36) * (Fraction(15) - Fraction(39, 22)) ** 2,
        "energy45": stationary_upper_dimension(13, 45) * (Fraction(15) - Fraction(39, 23)) ** 2,
    }
    assert all(value < q132 for value in size13_q132.values())

    # If a determinant above q132 used (14,1) on both axes, internally
    # isolated cases are impossible by the equal-color-sign orbit lemma.
    # A no-isolate energy-63 case is the friendship support.  At energy 72,
    # exact internal enumeration suffices.  At energy 81 the non-friendship
    # possibility is nine unit edges on 14 nonisolated vertices; components
    # with more than five edges cannot occur by the vertex-coverage count.
    e72_no_isolate = EXPECTED_MAXIMA[72] * (Fraction(15) - Fraction(51, 1) / LAMBDA_14[72])
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    e81_layer = q9_layer(connected_types(5, norm9), 9, required_isolates=0, order=14)
    assert len(e81_layer["profiles"]) == 18
    assert e81_layer["maximum_determinant"] == 20438941279518720
    e81_no_isolate = e81_layer["maximum_determinant"] * (
        Fraction(15) - Fraction(42, 1) / Fraction(261, 10)
    )
    assert e72_no_isolate < q132 and e81_no_isolate < q132

    return {
        "q126_envelope_numerator": q126.numerator,
        "q126_envelope_denominator": q126.denominator,
        "q123_color_14_1_splits": rows14,
        "q123_color_13_2_generic": generic_rows,
        "q123_color_13_2_tight": tight_rows,
        "q132_audit": {
            "size13_shape_over_q132": {key: float(value / q132) for key, value in size13_q132.items()},
            "energy72_no_isolate_over_q132": float(e72_no_isolate / q132),
            "energy81_no_isolate_profiles": len(e81_layer["profiles"]),
            "energy81_no_isolate_internal_maximum": e81_layer["maximum_determinant"],
            "energy81_no_isolate_over_q132": float(e81_no_isolate / q132),
        },
        "theorem": "Every Q=123 color case lies below the Q=132 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
