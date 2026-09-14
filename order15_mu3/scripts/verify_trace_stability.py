"""Exact rational certificate for the trace-G^2 counterexample energy bound."""

from __future__ import annotations

import json
from fractions import Fraction
from math import isqrt


ORDER = 15
BENCHMARK = 2**22 * 3**20 * 19


def sqrt_lower(value: Fraction, scale: int = 10**15) -> Fraction:
    """Return a certified rational lower bound for a positive square root."""

    numerator = isqrt(value.numerator * scale * scale // value.denominator)
    result = Fraction(numerator, scale)
    assert result * result <= value
    return result


def stationary_upper(energy: int, multiplicity: int) -> Fraction:
    """Upper-bound the two-level stationary product at fixed variance."""

    m = multiplicity
    t = sqrt_lower(Fraction(energy, ORDER * m * (ORDER - m)))
    high = Fraction(ORDER) + (ORDER - m) * t
    low = Fraction(ORDER) - m * t
    if low <= 0:
        return Fraction(0)
    return high**m * low ** (ORDER - m)


def main() -> None:
    # Q=sum_{i<j}|G_ij|^2 and S=tr((G-15I)^2)=2Q.
    q_ruled = 171
    energy = 2 * q_ruled
    candidates = {m: stationary_upper(energy, m) for m in range(1, ORDER)}
    maximizing_m = max(candidates, key=candidates.get)
    maximum = candidates[maximizing_m]
    assert maximum < BENCHMARK
    assert q_ruled % 3 == 0
    print(json.dumps({
        "first_excluded_total_off_diagonal_norm": q_ruled,
        "variance_energy": energy,
        "maximizing_stationary_multiplicity": maximizing_m,
        "certified_upper_numerator": maximum.numerator,
        "certified_upper_denominator": maximum.denominator,
        "certified_upper_over_benchmark": float(maximum / BENCHMARK),
        "counterexample_necessary_bound": "sum_{i<j} |G_ij|^2 <= 168",
    }, indent=2))


if __name__ == "__main__":
    main()
