from order15_mu3.scripts.verify_pre_q126_audit import certificate


def test_every_pre_q126_bound_is_below_q126() -> None:
    result = certificate()
    assert max(result["size13_exact_shape_over_q126"].values()) < 1
    assert max(result["size13_exact_shape_over_q132"].values()) < 1
    assert result["size14_no_isolate_cases"]["q114_e72_over_q126"] < 1
    assert result["size14_no_isolate_cases"]["q114_e72_over_q132"] < 1
