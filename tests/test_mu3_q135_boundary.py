from order15_mu3.scripts.verify_q135_boundary import certificate


def test_q135_is_below_q141_envelope() -> None:
    result = certificate()
    assert result["all_same_over_q141"] < 1
    assert max(result["size13_over_q141"].values()) < 1
