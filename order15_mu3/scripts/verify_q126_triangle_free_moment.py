"""Third-moment determinant certificate for triangle-free Q=126 supports."""

from __future__ import annotations

import json
from fractions import Fraction
from math import factorial


UPPER = Fraction(4 * 15**15, 7)
NO_K4_PURE_UPPER = Fraction(13 * 15**15, 20)


def certificate() -> dict[str, object]:
    # For E=G-15I, triangle-freeness gives tr(E^3)=0, while Q=126 gives
    # tr(E)=0 and tr(E^2)=252.  For every t>-1,
    #   log(1+t) <= t-t^2/2+t^3/3,
    # since the derivative of the difference is t^3/(1+t).
    # Summing at t=lambda(E)/15 yields det(G)<=15^15 exp(-14/25).
    # The exact Taylor lower bound below proves exp(14/25)>7/4.
    x = Fraction(14, 25)
    exponential_lower = sum((x**degree / factorial(degree) for degree in range(9)), Fraction(0))
    assert exponential_lower > Fraction(7, 4)
    # A K4-free graph with fourteen edges has at most eight triangles (the
    # sharp clique-density extremum is a subgraph of a complete tripartite
    # graph).  Each unit triangle contributes at most 162 to tr(E^3).
    exponent = Fraction(14, 25) - Fraction(8 * 162, 3 * 15**3)
    pure_lower = sum((exponent**degree / factorial(degree) for degree in range(7)), Fraction(0))
    assert pure_lower > Fraction(20, 13)
    return {
        "trace_E": 0,
        "trace_E2": 252,
        "trace_E3": 0,
        "exp_lower_numerator": exponential_lower.numerator,
        "exp_lower_denominator": exponential_lower.denominator,
        "determinant_upper_numerator": UPPER.numerator,
        "determinant_upper_denominator": UPPER.denominator,
        "no_K4_pure_determinant_upper_numerator": NO_K4_PURE_UPPER.numerator,
        "no_K4_pure_determinant_upper_denominator": NO_K4_PURE_UPPER.denominator,
        "theorem": "Triangle-free and pure K4-free Q=126 Grams satisfy the displayed cubic-log bounds.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
