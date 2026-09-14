from order15_mu3.scripts.verify_q87_boundary import certificate


def test_q87_boundary_is_below_q90_envelope() -> None:
    result = certificate()
    assert len(result["color_14_1_splits"]) == 4
    assert all(row["over_q90_envelope"] < 1 for row in result["color_14_1_splits"])
