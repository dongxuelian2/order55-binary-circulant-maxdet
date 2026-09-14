"""Aggregate exact branch bounds and apply the arithmetic determinant sieve."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_pre_q144_audit import certificate as pre_q144_certificate
    from order15_mu3.scripts.verify_q144_boundary import certificate as q144_certificate
    from order15_mu3.scripts.verify_q150_boundary import certificate as q150_certificate
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, stationary_upper
except ModuleNotFoundError:
    from verify_bounds import is_rational_eisenstein_norm
    from verify_pre_q144_audit import certificate as pre_q144_certificate
    from verify_q144_boundary import certificate as q144_certificate
    from verify_q150_boundary import certificate as q150_certificate
    from verify_trace_stability import BENCHMARK, stationary_upper


def certificate() -> dict[str, object]:
    earlier = pre_q144_certificate()
    boundary = q144_certificate()
    q150 = q150_certificate()
    earlier_upper = Fraction(
        earlier["certified_upper_numerator"], earlier["certified_upper_denominator"]
    )
    boundary_upper = Fraction(
        boundary["certified_upper_numerator"], boundary["certified_upper_denominator"]
    )
    q150_upper = Fraction(q150["certified_upper_numerator"], q150["certified_upper_denominator"])
    tail_values = {
        q: max(stationary_upper(2 * q, multiplicity) for multiplicity in range(1, 15))
        for q in range(153, 169, 3)
    }
    tail_q = max(tail_values, key=tail_values.get)
    assert tail_q == 153
    raw_upper = max(earlier_upper, boundary_upper, q150_upper, tail_values[tail_q])
    sources = {
        "Q<144 exact branch aggregate": earlier_upper,
        "Q=144 exact boundary": boundary_upper,
        "Q=150 exact boundary": q150_upper,
        "Q>=153 trace envelope": tail_values[tail_q],
    }
    source = max(sources, key=sources.get)

    # Every squared determinant is an integer divisible by 3^14 and a norm
    # from Z[omega].  Because all analytic certificates are strict, the
    # largest possible integer is ceil(raw_upper)-1.
    integer_ceiling_minus_one = (raw_upper.numerator - 1) // raw_upper.denominator
    divisor = 3**14
    candidate = integer_ceiling_minus_one // divisor * divisor
    skipped = 0
    while not is_rational_eisenstein_norm(candidate):
        candidate -= divisor
        skipped += 1
    assert candidate >= BENCHMARK
    return {
        "source_bounds": {
            name: {
                "numerator": value.numerator,
                "denominator": value.denominator,
                "decimal": float(value),
            }
            for name, value in sources.items()
        },
        "maximizing_source": source,
        "tail_maximizing_energy": tail_q,
        "raw_upper_numerator": raw_upper.numerator,
        "raw_upper_denominator": raw_upper.denominator,
        "integer_ceiling_minus_one": integer_ceiling_minus_one,
        "divisibility": "3^14",
        "non_norm_multiples_skipped": skipped,
        "rigorous_integer_upper": candidate,
        "benchmark": BENCHMARK,
        "squared_gap_factor": candidate / BENCHMARK,
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
