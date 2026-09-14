"""Aggregate exact branch bounds after the post-Q150 tail closure.

The high-energy shells Q=153,159,162,168 are now exact exclusions, not
analytic envelope candidates.  Therefore the global upper is controlled by
the surviving Q<=150 certificates only, followed by the same 3^14 divisibility
and rational Eisenstein-norm sieve used previously.
"""

from __future__ import annotations

import json
from fractions import Fraction
from math import sqrt

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_pre_q144_audit import certificate as pre_q144_certificate
    from order15_mu3.scripts.verify_q144_boundary import certificate as q144_certificate
    from order15_mu3.scripts.verify_q150_boundary import certificate as q150_certificate
    from order15_mu3.scripts.verify_tail_final_closure import certificate as tail_closure_certificate
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:
    from verify_bounds import is_rational_eisenstein_norm
    from verify_pre_q144_audit import certificate as pre_q144_certificate
    from verify_q144_boundary import certificate as q144_certificate
    from verify_q150_boundary import certificate as q150_certificate
    from verify_tail_final_closure import certificate as tail_closure_certificate
    from verify_trace_stability import BENCHMARK


HADAMARD_SQUARED = 15**15


def normalized_percent(squared_determinant: int) -> float:
    """Return 100*|det(H)|/15^(15/2) from a squared determinant."""

    return 100.0 * sqrt(squared_determinant / HADAMARD_SQUARED)


def certificate() -> dict[str, object]:
    earlier = pre_q144_certificate()
    boundary = q144_certificate()
    q150 = q150_certificate()
    tail = tail_closure_certificate()
    assert tail["closed_tail_shells"] == [153, 159, 162, 168]
    assert tail["new_counterexample_energy_bound"] == "Q <= 150"

    earlier_upper = Fraction(
        earlier["certified_upper_numerator"], earlier["certified_upper_denominator"]
    )
    boundary_upper = Fraction(
        boundary["certified_upper_numerator"], boundary["certified_upper_denominator"]
    )
    q150_upper = Fraction(
        q150["certified_upper_numerator"], q150["certified_upper_denominator"]
    )

    # The old aggregate also included the Q=153 trace envelope.  That is no
    # longer a valid controlling source because every genuine Q>150 shell is
    # now structurally excluded below the benchmark.
    sources = {
        "Q<144 exact branch aggregate": earlier_upper,
        "Q=144 exact boundary": boundary_upper,
        "Q=150 exact boundary": q150_upper,
    }
    source = max(sources, key=sources.get)
    raw_upper = sources[source]

    # Every squared determinant is an integer divisible by 3^14 and a norm
    # from Z[omega].  All analytic source bounds are strict, so the largest
    # possible integer is ceil(raw_upper)-1 before the arithmetic sieve.
    integer_ceiling_minus_one = (raw_upper.numerator - 1) // raw_upper.denominator
    divisor = 3**14
    candidate = integer_ceiling_minus_one // divisor * divisor
    skipped = 0
    while not is_rational_eisenstein_norm(candidate):
        candidate -= divisor
        skipped += 1

    assert candidate >= BENCHMARK
    assert candidate < 287710229992756239
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
        "closed_post_q150_shells": tail["closed_tail_shells"],
        "counterexample_energy_bound": tail["new_counterexample_energy_bound"],
        "raw_upper_numerator": raw_upper.numerator,
        "raw_upper_denominator": raw_upper.denominator,
        "integer_ceiling_minus_one": integer_ceiling_minus_one,
        "divisibility": "3^14",
        "non_norm_multiples_skipped": skipped,
        "rigorous_integer_upper": candidate,
        "benchmark": BENCHMARK,
        "squared_gap_factor": candidate / BENCHMARK,
        "normalized_lower_percent": normalized_percent(BENCHMARK),
        "normalized_upper_percent": normalized_percent(candidate),
        "normalized_gap_percentage_points": normalized_percent(candidate) - normalized_percent(BENCHMARK),
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
