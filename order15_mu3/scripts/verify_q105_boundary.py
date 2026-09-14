"""Exact Schur certificate placing Q=99 and Q=105 below the Q=108 envelope."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_15_energy99 import Q99_ALL_SAME_UPPER
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_15_energy99 import Q99_ALL_SAME_UPPER
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_trace_stability import stationary_upper


LAMBDA_14 = {
    18: Fraction(104, 5),
    27: Fraction(221, 10),
    36: Fraction(116, 5),
    45: Fraction(121, 5),
    54: Fraction(251, 10),
    # At energy 63 there are at most seven edges, hence clique number <=4.
    # Cauchy--Schwarz plus Motzkin--Straus gives rho(E)^2<=2e(1-1/4)=189/2.
    63: Fraction(99, 4),
}
LAMBDA_13 = {9: Fraction(96, 5), 18: Fraction(21), 27: Fraction(221, 10)}


def certificate() -> dict[str, object]:
    q108_envelope = max(stationary_upper(216, multiplicity) for multiplicity in range(1, 15))
    rows_14 = []
    for internal_energy, lambda_upper in LAMBDA_14.items():
        cross_energy = 105 - internal_energy
        if internal_energy == 63:
            assert (lambda_upper - 15) ** 2 > Fraction(189, 2)
        else:
            assert (lambda_upper - 15) ** 2 > Fraction(2 * internal_energy * 13, 14)
        upper = stationary_upper_dimension(14, internal_energy) * (
            15 - Fraction(cross_energy, 1) / lambda_upper
        )
        assert upper < q108_envelope
        rows_14.append({
            "internal_energy": internal_energy,
            "cross_energy": cross_energy,
            "lambda_upper": [lambda_upper.numerator, lambda_upper.denominator],
            "over_q108_envelope": float(upper / q108_envelope),
        })

    rows_13 = []
    for internal_energy, lambda_upper in LAMBDA_13.items():
        assert (lambda_upper - 15) ** 2 > Fraction(2 * internal_energy * 12, 13)
        upper = stationary_upper_dimension(13, internal_energy) * (
            15 - Fraction(78, 2) / lambda_upper
        ) ** 2
        assert upper < q108_envelope
        rows_13.append({
            "internal_energy": internal_energy,
            "cross_energy_lower_bound": 78,
            "lambda_upper": [lambda_upper.numerator, lambda_upper.denominator],
            "over_q108_envelope": float(upper / q108_envelope),
        })

    q99_size13_upper = rows_13[0]["over_q108_envelope"]
    assert Q99_ALL_SAME_UPPER < q108_envelope
    return {
        "q99_all_same_upper": Q99_ALL_SAME_UPPER,
        "q99_size13_over_q108_envelope": q99_size13_upper,
        "q105_color_14_1_splits": rows_14,
        "q105_color_13_2_splits": rows_13,
        "q108_envelope_numerator": q108_envelope.numerator,
        "q108_envelope_denominator": q108_envelope.denominator,
        "theorem": "Every Q=99 and Q=105 counterexample lies below the Q=108 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
