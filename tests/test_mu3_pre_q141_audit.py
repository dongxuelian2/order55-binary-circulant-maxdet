from order15_mu3.scripts.verify_pre_q141_audit import certificate


def test_q126_k5_tradeoff_is_below_q141() -> None:
    result = certificate()
    assert result["no_k5_over_q141"] < 1
    assert max(result["weighted_k5_over_q141"].values()) < 1
    assert max(result["pure_k5_tradeoff_over_q141"].values()) < 1
