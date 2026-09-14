from order15_mu3.scripts.verify_q126_boundary import certificate


def test_q126_is_below_q132_envelope() -> None:
    result = certificate()
    assert result["all_same_over_q132"] < 1
    assert len(result["color_13_1_1_generic"]) == 9
    assert len(result["color_13_1_1_tight"]) == 4
    assert all(row["over_q132"] < 1 for row in result["color_13_1_1_generic"])
    assert all(row["over_q132"] < 1 for row in result["color_13_1_1_tight"])
