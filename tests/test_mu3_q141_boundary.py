from order15_mu3.scripts.verify_q141_boundary import certificate


def test_q141_is_below_q144_envelope() -> None:
    result = certificate()
    assert result["size13_max_over_q144"] < 1
    assert max(result["size14_over_q144"].values()) < 1
