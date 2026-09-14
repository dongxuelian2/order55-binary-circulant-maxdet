from fractions import Fraction

from order15_mu3.scripts.verify_q126_triangle_free_moment import NO_K4_PURE_UPPER, UPPER, certificate


def test_q126_triangle_free_third_moment_bound() -> None:
    result = certificate()
    assert Fraction(
        result["determinant_upper_numerator"], result["determinant_upper_denominator"]
    ) == UPPER
    assert UPPER < 278_000_000_000_000_000
    assert NO_K4_PURE_UPPER < 285_000_000_000_000_000
