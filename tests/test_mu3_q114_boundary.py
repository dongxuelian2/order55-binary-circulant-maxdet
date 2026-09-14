from order15_mu3.scripts.verify_q114_boundary import certificate


def test_q108_and_q114_are_below_q117_envelope() -> None:
    result = certificate()
    for key in ("q108_color_13_1_1_cases", "q114_color_14_1_cases", "q114_color_13_2_cases"):
        assert result[key]
        assert all(row["over_q117"] < 1 for row in result[key])
