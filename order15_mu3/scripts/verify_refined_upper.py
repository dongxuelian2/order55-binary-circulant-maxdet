"""Exact certificate for the refined global order-15 mu3 upper bound."""

from __future__ import annotations

import json

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK, stationary_upper
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_bounds import is_rational_eisenstein_norm
    from verify_trace_stability import BENCHMARK, stationary_upper


Q78_SLICE_MAXIMA = {
    "internal18_cross60": 298188223593750000,
    "internal27_cross51": 307903629375000000,
    "internal36_cross42": 308697676962890625,
}
Q87_SIZE13_MAXIMUM = 289967495625000000
Q90_BOUND = 289967495625000000
Q99_ALL_SAME_BOUND = 282438239494864896
REFINED_UPPER = 293695860429200709


def certificate() -> dict[str, object]:
    # Q=60 and 69 are excluded by the twin-leaf certificate.  Q=78 is
    # excluded by the equal-color-sign isolated-orbit obstruction.  At Q=87, the
    # (13,2,0) exact maximum and four (14,1,0) Schur bounds lie below the Q=90
    # envelope.  The exact Q=90 support certificate lies below the Q=99
    # envelope, and every Q=96 color split does too.  The Q=99 all-same
    # support certificate and its size-13 Schur case lie below the Q=108
    # envelope; every Q=105 split does too.  The all-same Q=108 class is
    # excluded above the record, while its size-13 class and all Q=114 splits
    # lie below the Q=117 envelope (verify_q114_boundary.py).  Exact internal
    # support refinements, the same color-orbit obstruction, the friendship
    # theorem, and exact support/Schur certificates at Q=117 and Q=123 then
    # put every energy below 126 under the Q=126 envelope
    # (verify_pre_q126_audit.py).  The later support/Schur certificates close
    # Q=126,132,135,141 below the Q=144 envelope; in particular the Q=126
    # no-K5 case is completed by the K4-component energy/isolate tradeoff in
    # verify_pre_q144_audit.py.
    # Energies congruent to 3 mod 9 are included harmlessly although color
    # congruence makes those slices empty.
    trace_values = {
        q: max(stationary_upper(2 * q, multiplicity) for multiplicity in range(1, 15))
        for q in range(144, 169, 3)
    }
    maximum_q = max(trace_values, key=trace_values.get)
    continuous = trace_values[maximum_q]
    assert maximum_q == 144
    assert Q87_SIZE13_MAXIMUM < continuous
    assert Q90_BOUND < continuous
    assert Q99_ALL_SAME_BOUND < continuous
    integer_floor = continuous.numerator // continuous.denominator

    divisor = 3**14
    candidate = integer_floor // divisor * divisor
    skipped = 0
    while not is_rational_eisenstein_norm(candidate):
        candidate -= divisor
        skipped += 1
    assert skipped == 4 and candidate == REFINED_UPPER
    return {
        "q78_arithmetic_slice_maxima": Q78_SLICE_MAXIMA,
        "q87_size13_exact_maximum": Q87_SIZE13_MAXIMUM,
        "q90_certified_upper": Q90_BOUND,
        "q99_all_same_certified_upper": Q99_ALL_SAME_BOUND,
        "trace_range_checked": [144, 168],
        "maximizing_energy": maximum_q,
        "continuous_upper_numerator": continuous.numerator,
        "continuous_upper_denominator": continuous.denominator,
        "continuous_integer_floor": integer_floor,
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
