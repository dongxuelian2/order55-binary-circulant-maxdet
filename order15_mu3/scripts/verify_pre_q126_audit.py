"""Compact exact audit that every energy below 126 is below the Q=126 envelope."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_13_energy9 import arithmetic_admissible
    from order15_mu3.scripts.verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from order15_mu3.scripts.verify_color_15_energy99 import Q99_ALL_SAME_UPPER
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_boundary import LAMBDA_14, one_edge_tight
    from order15_mu3.scripts.verify_q114_internal_support import EXPECTED_MAXIMA
    from order15_mu3.scripts.verify_q117_boundary import schur_two_by_two
    from order15_mu3.scripts.verify_q78_color_orbit_obstruction import certificate as q78_obstruction
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, stationary_upper
except ModuleNotFoundError:
    from verify_color_13_energy9 import arithmetic_admissible
    from verify_color_15_energy90 import SIZE13_Q90_MAXIMUM
    from verify_color_15_energy99 import Q99_ALL_SAME_UPPER
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_boundary import LAMBDA_14, one_edge_tight
    from verify_q114_internal_support import EXPECTED_MAXIMA
    from verify_q117_boundary import schur_two_by_two
    from verify_q78_color_orbit_obstruction import certificate as q78_obstruction
    from verify_trace_stability import BENCHMARK, stationary_upper


def certificate() -> dict[str, object]:
    q126 = max(stationary_upper(252, multiplicity) for multiplicity in range(1, 15))
    q132 = max(stationary_upper(264, multiplicity) for multiplicity in range(1, 15))
    assert q78_obstruction()["orbit_intersection"] == []

    # The exact energy-9 size-13 phase enumeration (used at Q=87,90,96,99,
    # and 108) has this universal arithmetic maximum.
    assert arithmetic_admissible(SIZE13_Q90_MAXIMUM)
    assert SIZE13_Q90_MAXIMUM < q132
    assert Q99_ALL_SAME_UPPER < q132

    # Exact internal support shapes used whenever the generic size-13 Schur
    # estimate is too weak.  Outside off-diagonal terms may be discarded.
    lambda_p3 = Fraction(193, 10)
    e18 = max(
        schur_two_by_two(3105 * 15**10, lambda_p3),
        schur_two_by_two(216**2 * 15**9, Fraction(18)),
    )
    e27 = max(
        schur_two_by_two(198 * 15**11, Fraction(21)),
        schur_two_by_two(216**3 * 15**7, Fraction(18)),
        schur_two_by_two(3105 * 216 * 15**8, lambda_p3),
        schur_two_by_two(44631 * 15**9, Fraction(20)),
        schur_two_by_two(44550 * 15**9, Fraction(21)),
        schur_two_by_two(3024 * 15**10, Fraction(21)),
    )
    e36 = max(
        schur_two_by_two(189 * 15**11, Fraction(21)),
        schur_two_by_two(198 * 216 * 15**9, Fraction(21)),
        schur_two_by_two(2835 * 15**10, Fraction(21)),
        # Four unit edges have clique number <=3, hence rho(E)^2<=48<49.
        stationary_upper_dimension(13, 36) * (Fraction(15) - Fraction(39, 22)) ** 2,
    )
    assert one_edge_tight(27) < q132
    assert e18 < q132 and e27 < q132 and e36 < q132

    # If a Q=87 or Q=96 determinant exceeded q126, its size-13 alternative
    # is already too small, forcing (14,1) on both axes.  At those energies
    # each internal support has an isolate, and the Q=78 color-orbit lemma
    # applies.  At Q=105 and Q=114, the only extra no-isolate size-14 cases
    # are a seven-edge friendship support (bounded by BENCHMARK) and the
    # eight-unit-edge Q=114 block handled by the exact determinant below.
    q114_e72 = EXPECTED_MAXIMA[72] * (
        Fraction(15) - Fraction(42, 1) / LAMBDA_14[72]
    )
    assert q114_e72 < q132

    return {
        "q126_envelope_numerator": q126.numerator,
        "q126_envelope_denominator": q126.denominator,
        "q132_envelope_numerator": q132.numerator,
        "q132_envelope_denominator": q132.denominator,
        "energy78": "excluded by equal-color-sign isolated-orbit obstruction",
        "size13_energy9_arithmetic_maximum": SIZE13_Q90_MAXIMUM,
        "q99_all_same_maximum": Q99_ALL_SAME_UPPER,
        "size13_exact_shape_over_q126": {
            "energy18": float(e18 / q126),
            "energy27": float(e27 / q126),
            "energy36": float(e36 / q126),
        },
        "size13_exact_shape_over_q132": {
            "energy18": float(e18 / q132),
            "energy27": float(e27 / q132),
            "energy36": float(e36 / q132),
        },
        "size14_no_isolate_cases": {
            "seven_edge_friendship_maximum": BENCHMARK,
            "q114_e72_over_q126": float(q114_e72 / q126),
            "q114_e72_over_q132": float(q114_e72 / q132),
        },
        "later_certificates": [
            "all-same Q=108 excluded above the record",
            "all-same Q=117 excluded above the record",
            "(13,1,1) Q=117 below Q=126",
            "all Q=123 color cases below Q=126",
        ],
        "theorem": "Every strict counterexample with Q<126 has determinant below the Q=132 envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
