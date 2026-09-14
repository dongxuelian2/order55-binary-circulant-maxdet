"""Exact Schur certificate placing every Q=87 case below the Q=90 envelope."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, stationary_upper
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_trace_stability import BENCHMARK, stationary_upper


SIZE13_BOUNDARY_MAXIMUM = 289967495625000000
LAMBDA_BOUNDS = {
    18: Fraction(104, 5),
    27: Fraction(221, 10),
    36: Fraction(116, 5),
    45: Fraction(121, 5),
}


def certificate() -> dict[str, object]:
    q90_envelope = max(stationary_upper(180, multiplicity) for multiplicity in range(1, 15))
    rows = []
    for internal_energy, lambda_upper in LAMBDA_BOUNDS.items():
        cross_energy = 87 - internal_energy
        # lambda_max <= 15 + sqrt(2 e * 13/14).
        assert (lambda_upper - 15) ** 2 > Fraction(2 * internal_energy * 13, 14)
        upper = stationary_upper_dimension(14, internal_energy) * (15 - Fraction(cross_energy, 1) / lambda_upper)
        assert upper < q90_envelope
        rows.append({
            "internal_energy": internal_energy,
            "cross_energy": cross_energy,
            "lambda_upper": [lambda_upper.numerator, lambda_upper.denominator],
            "schur_upper_numerator": upper.numerator,
            "schur_upper_denominator": upper.denominator,
            "over_q90_envelope": float(upper / q90_envelope),
        })
    assert SIZE13_BOUNDARY_MAXIMUM < q90_envelope
    return {
        "total_energy": 87,
        "color_14_1_splits": rows,
        "color_13_2_exact_arithmetic_maximum": SIZE13_BOUNDARY_MAXIMUM,
        "q90_envelope_numerator": q90_envelope.numerator,
        "q90_envelope_denominator": q90_envelope.denominator,
        "q90_envelope_over_record": float(q90_envelope / BENCHMARK),
        "theorem": "Every Q=87 counterexample lies below the Q=90 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
