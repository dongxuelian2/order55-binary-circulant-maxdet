"""Exact Schur certificate placing every Q=96 case below the Q=99 envelope."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_trace_stability import stationary_upper


LAMBDA_14 = {
    18: Fraction(104, 5),
    27: Fraction(221, 10),
    36: Fraction(116, 5),
    45: Fraction(121, 5),
    54: Fraction(251, 10),
}
LAMBDA_13 = {9: Fraction(96, 5), 18: Fraction(21)}


def certificate() -> dict[str, object]:
    q99_envelope = max(stationary_upper(198, multiplicity) for multiplicity in range(1, 15))
    rows_14 = []
    for internal_energy, lambda_upper in LAMBDA_14.items():
        cross_energy = 96 - internal_energy
        assert (lambda_upper - 15) ** 2 > Fraction(2 * internal_energy * 13, 14)
        upper = stationary_upper_dimension(14, internal_energy) * (
            15 - Fraction(cross_energy, 1) / lambda_upper
        )
        assert upper < q99_envelope
        rows_14.append({
            "internal_energy": internal_energy,
            "cross_energy": cross_energy,
            "lambda_upper": [lambda_upper.numerator, lambda_upper.denominator],
            "upper_numerator": upper.numerator,
            "upper_denominator": upper.denominator,
            "over_q99_envelope": float(upper / q99_envelope),
        })

    # In (13,2,0), cross energy is at least 78.  If the large block has
    # energy 9 the outside pair may carry the remaining 9; discarding that
    # off-diagonal entry only enlarges the 2x2 Schur determinant.  Thus c=78
    # safely dominates both subcases.
    rows_13 = []
    for internal_energy, lambda_upper in LAMBDA_13.items():
        cross_energy = 78
        assert (lambda_upper - 15) ** 2 > Fraction(2 * internal_energy * 12, 13)
        upper = stationary_upper_dimension(13, internal_energy) * (
            15 - Fraction(cross_energy, 2) / lambda_upper
        ) ** 2
        assert upper < q99_envelope
        rows_13.append({
            "internal_energy": internal_energy,
            "cross_energy_lower_bound": cross_energy,
            "lambda_upper": [lambda_upper.numerator, lambda_upper.denominator],
            "upper_numerator": upper.numerator,
            "upper_denominator": upper.denominator,
            "over_q99_envelope": float(upper / q99_envelope),
        })
    return {
        "total_energy": 96,
        "color_14_1_splits": rows_14,
        "color_13_2_splits": rows_13,
        "q99_envelope_numerator": q99_envelope.numerator,
        "q99_envelope_denominator": q99_envelope.denominator,
        "theorem": "Every Q=96 counterexample lies below the Q=99 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
