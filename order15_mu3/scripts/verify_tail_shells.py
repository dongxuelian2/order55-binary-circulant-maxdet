"""Shellwise trace/arithmetic certificates for the genuine high-energy tail.

The color congruence leaves exactly Q in {153,159,162,168} after the Q=150
boundary certificate and before trace stability excludes Q>=171.  This module
certifies each remaining shell separately rather than padding the tail with the
impossible Q=156 and Q=165 values.
"""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, stationary_upper
except ModuleNotFoundError:
    from verify_bounds import is_rational_eisenstein_norm
    from verify_trace_stability import BENCHMARK, stationary_upper


TAIL_SHELLS = (153, 159, 162, 168)
DIVISOR = 3**14


def shell_certificate(q: int) -> dict[str, object]:
    assert q in TAIL_SHELLS
    candidates = {m: stationary_upper(2 * q, m) for m in range(1, 15)}
    maximizing_m = max(candidates, key=candidates.get)
    envelope = candidates[maximizing_m]

    # Analytic bounds are strict, so the largest possible integer square is
    # ceil(envelope)-1.  Squared determinants are divisible by 3^14 and are
    # rational Eisenstein norms.
    integer_ceiling_minus_one = (envelope.numerator - 1) // envelope.denominator
    arithmetic_candidate = integer_ceiling_minus_one // DIVISOR * DIVISOR
    skipped = 0
    while arithmetic_candidate > 0 and not is_rational_eisenstein_norm(arithmetic_candidate):
        arithmetic_candidate -= DIVISOR
        skipped += 1

    return {
        "Q": q,
        "variance_energy": 2 * q,
        "maximizing_stationary_multiplicity": maximizing_m,
        "envelope_numerator": envelope.numerator,
        "envelope_denominator": envelope.denominator,
        "envelope_over_benchmark": float(envelope / BENCHMARK),
        "integer_ceiling_minus_one": integer_ceiling_minus_one,
        "non_norm_multiples_skipped": skipped,
        "largest_arithmetic_survivor": arithmetic_candidate,
        "arithmetic_survivor_over_benchmark": arithmetic_candidate / BENCHMARK,
        "record_excluded_by_shell_bound": arithmetic_candidate < BENCHMARK,
    }


def certificate() -> dict[str, object]:
    rows = [shell_certificate(q) for q in TAIL_SHELLS]
    assert all(row["maximizing_stationary_multiplicity"] == 1 for row in rows)
    envelopes = [Fraction(row["envelope_numerator"], row["envelope_denominator"]) for row in rows]
    assert envelopes == sorted(envelopes, reverse=True)
    assert rows[0]["Q"] == 153
    assert rows[-1]["Q"] == 168
    assert rows[-1]["envelope_over_benchmark"] > 1
    assert max(stationary_upper(342, m) for m in range(1, 15)) < BENCHMARK
    return {
        "genuine_tail_shells": list(TAIL_SHELLS),
        "shells": rows,
        "tail_maximizing_Q": 153,
        "tail_maximizing_stationary_multiplicity": 1,
        "first_trace_excluded_Q": 171,
        "theorem": (
            "The only genuine post-Q150 shells before trace stability are "
            "Q=153,159,162,168; their stationary envelopes decrease strictly, "
            "are all maximized at multiplicity one, and Q>=171 is below the record."
        ),
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
