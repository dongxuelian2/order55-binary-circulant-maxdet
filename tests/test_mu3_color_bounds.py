from fractions import Fraction

from order15_mu3.scripts.verify_color_partition_bounds import schur_upper
from order15_mu3.scripts.verify_trace_stability import BENCHMARK


def test_decisive_color_schur_cases() -> None:
    assert schur_upper(10, 9, Fraction(18), 150) < BENCHMARK
    assert schur_upper(11, 27, Fraction(23), 132) < BENCHMARK
    assert schur_upper(12, 27, Fraction(21), 108) < BENCHMARK
    assert schur_upper(12, 54, Fraction(25), 108) < BENCHMARK
