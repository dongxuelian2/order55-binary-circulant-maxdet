from order15_mu3.scripts.verify_q144_boundary import certificate


def test_q144_is_below_q150_envelope() -> None:
    result = certificate()
    assert result["all_same_over_q150"] < 1
    assert result["size13_max_over_q150"] < 1
