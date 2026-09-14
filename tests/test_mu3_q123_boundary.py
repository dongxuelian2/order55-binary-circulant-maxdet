from order15_mu3.scripts.verify_q123_boundary import certificate


def test_q123_is_below_q126_envelope() -> None:
    result = certificate()
    assert len(result["q123_color_14_1_splits"]) == 9
    assert len(result["q123_color_13_2_generic"]) == 9
    assert len(result["q123_color_13_2_tight"]) == 3
    for key in ("q123_color_14_1_splits", "q123_color_13_2_generic", "q123_color_13_2_tight"):
        assert all(row["over_q126"] < 1 for row in result[key])
    audit = result["q132_audit"]
    assert max(audit["size13_shape_over_q132"].values()) < 1
    assert audit["energy72_no_isolate_over_q132"] < 1
    assert audit["energy81_no_isolate_over_q132"] < 1
