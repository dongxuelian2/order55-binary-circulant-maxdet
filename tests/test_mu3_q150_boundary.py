from order15_mu3.scripts.verify_q150_boundary import certificate


def test_q150_is_below_q153_envelope() -> None:
    result = certificate()
    assert result["size13_max_over_q153"] < 1
    assert max(result["size14_over_q153"].values()) < 1
