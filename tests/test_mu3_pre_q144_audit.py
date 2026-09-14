from order15_mu3.scripts.verify_pre_q144_audit import certificate


def test_every_earlier_slice_is_below_q144() -> None:
    result = certificate()
    assert max(result["size14_over_q144"].values()) < 1
    assert result["q126_no_k4_over_q144"] < 1
    assert result["q126_k4_max_over_q144"] < 1
    assert max(result["q126_weighted_k5_over_q144"].values()) < 1
    assert max(result["q126_pure_k5_over_q144"].values()) < 1
    assert result["q135_all_same_over_q144"] < 1
