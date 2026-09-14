from order15_mu3.scripts.verify_pre_q135_audit import certificate


def test_every_pre_q135_bound_is_below_q135() -> None:
    result = certificate()
    assert max(row["over_q135"] for row in result["size13_energy36"].values()) < 1
    assert max(result["size13_high_internal_over_q135"].values()) < 1
    assert result["q126_all_same_over_q135"] < 1
    assert max(result["size14_no_isolate_over_q135"].values()) < 1
