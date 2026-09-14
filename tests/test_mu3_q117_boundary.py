from order15_mu3.scripts.verify_q117_boundary import certificate


def test_q117_size13_is_below_q126_envelope() -> None:
    result = certificate()
    assert len(result["generic_allocations"]) == 6
    assert len(result["tight_allocations"]) == 4
    for key in ("generic_allocations", "tight_allocations"):
        assert all(row["over_q126"] < 1 for row in result[key])
    for key in ("energy18_shape_bounds", "energy27_shape_bounds", "energy36_shape_bounds"):
        assert result[key]
        assert max(result[key].values()) < 1
