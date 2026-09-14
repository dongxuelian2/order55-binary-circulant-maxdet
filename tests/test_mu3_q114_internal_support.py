from order15_mu3.scripts.verify_q114_internal_support import EXPECTED_MAXIMA, certificate


def test_q114_internal_support_maxima() -> None:
    result = certificate()
    assert {row["internal_energy"]: row["maximum_determinant"] for row in result["energy_rows"]} == EXPECTED_MAXIMA
