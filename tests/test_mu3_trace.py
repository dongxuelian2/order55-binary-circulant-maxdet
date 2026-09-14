from fractions import Fraction

from order15_mu3.scripts.verify_trace_stability import BENCHMARK, sqrt_lower, stationary_upper


def test_certified_square_root_lower_bound() -> None:
    value = Fraction(57, 35)
    lower = sqrt_lower(value)
    assert lower * lower <= value
    assert (lower + Fraction(1, 10**15)) ** 2 > value


def test_energy_171_is_below_record_at_every_stationary_type() -> None:
    assert max(stationary_upper(342, m) for m in range(1, 15)) < BENCHMARK
