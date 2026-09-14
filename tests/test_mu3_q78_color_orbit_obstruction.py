from order15_mu3.scripts.verify_q78_color_orbit_obstruction import certificate


def test_q78_color_orbit_obstruction() -> None:
    result = certificate()
    assert result["forced_row_and_column_partition"] == [14, 1, 0]
    assert result["minimum_internal_isolates_per_side"] >= 1
    assert result["orbit_intersection"] == []
