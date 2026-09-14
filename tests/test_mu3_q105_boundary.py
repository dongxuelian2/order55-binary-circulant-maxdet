from order15_mu3.scripts.verify_q105_boundary import certificate


def test_q99_and_q105_are_below_q108_envelope() -> None:
    result = certificate()
    assert result["q99_size13_over_q108_envelope"] < 1
    assert all(row["over_q108_envelope"] < 1 for row in result["q105_color_14_1_splits"])
    assert all(row["over_q108_envelope"] < 1 for row in result["q105_color_13_2_splits"])
