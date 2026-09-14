from order15_mu3.scripts.verify_q132_boundary import certificate


def test_q132_is_below_q135_envelope() -> None:
    result = certificate()
    assert max(result["size13_energy18_over_q135"].values()) < 1
    assert max(result["size13_energy27_over_q135"].values()) < 1
    assert all(row["over_q135"] < 1 for row in result["size13_generic"])
    assert max(result["size14_no_isolate_over_q135"].values()) < 1
