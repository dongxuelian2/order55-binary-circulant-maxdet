from order15_mu3.scripts.verify_q96_boundary import certificate


def test_q96_boundary_is_below_q99_envelope() -> None:
    result = certificate()
    assert len(result["color_14_1_splits"]) == 5
    assert len(result["color_13_2_splits"]) == 2
    assert all(row["over_q99_envelope"] < 1 for row in result["color_14_1_splits"])
    assert all(row["over_q99_envelope"] < 1 for row in result["color_13_2_splits"])
