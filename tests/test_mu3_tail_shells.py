from fractions import Fraction

from order15_mu3.scripts.verify_tail_shells import TAIL_SHELLS, certificate


def test_genuine_tail_shells_are_certified_individually() -> None:
    cert = certificate()
    assert tuple(cert["genuine_tail_shells"]) == TAIL_SHELLS == (153, 159, 162, 168)
    rows = cert["shells"]
    assert [row["Q"] for row in rows] == [153, 159, 162, 168]
    assert all(row["maximizing_stationary_multiplicity"] == 1 for row in rows)
    envelopes = [Fraction(row["envelope_numerator"], row["envelope_denominator"]) for row in rows]
    assert envelopes[0] > envelopes[1] > envelopes[2] > envelopes[3]
    assert cert["tail_maximizing_Q"] == 153
    assert cert["first_trace_excluded_Q"] == 171
