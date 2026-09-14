from order15_mu3.scripts.verify_refined_upper import REFINED_UPPER, certificate


def test_refined_upper_certificate() -> None:
    result = certificate()
    assert result["rigorous_integer_upper"] == REFINED_UPPER
    assert result["squared_gap_factor"] < 1.057
