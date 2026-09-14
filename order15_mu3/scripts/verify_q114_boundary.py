"""Exact Schur certificate placing Q=108 and Q=114 below Q=117."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_trace_stability import sqrt_lower, stationary_upper
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_trace_stability import sqrt_lower, stationary_upper


LAMBDA_13 = {9: Fraction(96, 5), 18: Fraction(21), 27: Fraction(221, 10), 36: Fraction(47, 2)}
LAMBDA_14 = {
    18: Fraction(104, 5), 27: Fraction(221, 10), 36: Fraction(116, 5),
    45: Fraction(121, 5), 54: Fraction(251, 10), 63: Fraction(99, 4),
    # Eight edges have clique number at most four, so rho(E)^2<=108.
    72: Fraction(127, 5),
}


def generic_size13(internal_energy: int, cross_energy: int) -> Fraction:
    bound = LAMBDA_13[internal_energy]
    assert (bound - 15) ** 2 > Fraction(2 * internal_energy * 12, 13)
    return stationary_upper_dimension(13, internal_energy) * (
        15 - Fraction(cross_energy, 2) / bound
    ) ** 2


def one_edge_tight(outside_norm: int) -> Fraction:
    """Energy-9 large block, cross energy 78, and mandatory outside edge."""

    block_determinant = 216 * 15**11
    trace_upper = Fraction(30) - Fraction(78, 18)
    correction_upper = Fraction(39, 12)
    residual_lower = sqrt_lower(Fraction(outside_norm)) - correction_upper
    assert residual_lower > 0
    return block_determinant * ((trace_upper / 2) ** 2 - residual_lower**2)


def certificate() -> dict[str, object]:
    q117 = max(stationary_upper(234, multiplicity) for multiplicity in range(1, 15))

    # Q=108, partition (13,1,1).  Only e=9,c=78,r=21 and
    # e=18,c=78,r=12 evade the generic Q117 comparison.
    q108_rows = []
    tight9 = one_edge_tight(21)
    assert tight9 < q117
    q108_rows.append({"internal": 9, "cross": 78, "outside_norm": 21, "over_q117": float(tight9 / q117)})
    # At internal energy 18 the support is P3 or 2K2.  Their exact block
    # determinants and lambda maxima give these trace-only bounds.
    adjacent_lambda = Fraction(193, 10)  # > 15+3 sqrt(2)
    assert (adjacent_lambda - 15) ** 2 > 18
    exact18 = (
        3105 * 15**10 * (Fraction(30) - Fraction(78, 1) / adjacent_lambda) ** 2 / 4,
        216**2 * 15**9 * (Fraction(30) - Fraction(78, 18)) ** 2 / 4,
    )
    assert all(value < q117 for value in exact18)
    q108_rows.append({"internal": 18, "cross": 78, "outside_norm": 12,
                      "over_q117": max(float(value / q117) for value in exact18)})
    for internal, cross in ((9, 87), (9, 96), (18, 87), (27, 78)):
        upper = generic_size13(internal, cross)
        assert upper < q117
        q108_rows.append({"internal": internal, "cross": cross, "over_q117": float(upper / q117)})

    rows14 = []
    for internal, bound in LAMBDA_14.items():
        cross = 114 - internal
        if internal == 63:
            assert (bound - 15) ** 2 > Fraction(189, 2)
        elif internal == 72:
            assert (bound - 15) ** 2 > 108
        else:
            assert (bound - 15) ** 2 > Fraction(2 * internal * 13, 14)
        upper = stationary_upper_dimension(14, internal) * (15 - Fraction(cross, 1) / bound)
        assert upper < q117
        rows14.append({"internal": internal, "cross": cross, "over_q117": float(upper / q117)})

    rows13 = []
    tight114 = one_edge_tight(27)
    assert tight114 < q117
    rows13.append({"internal": 9, "cross": 78, "outside_norm": 27, "over_q117": float(tight114 / q117)})
    for internal, cross in ((9, 96), (18, 87), (27, 78), (9, 105), (18, 96), (27, 87), (36, 78)):
        upper = generic_size13(internal, cross)
        assert upper < q117
        rows13.append({"internal": internal, "cross": cross, "over_q117": float(upper / q117)})
    return {
        "q108_color_13_1_1_cases": q108_rows,
        "q114_color_14_1_cases": rows14,
        "q114_color_13_2_cases": rows13,
        "q117_envelope_numerator": q117.numerator,
        "q117_envelope_denominator": q117.denominator,
        "theorem": "Every Q=108 and Q=114 counterexample lies below the Q=117 envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
