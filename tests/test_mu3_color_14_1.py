from order15_mu3.scripts.verify_color_14_1_low_energy import (
    one_internal_edge_candidates,
    two_internal_edge_candidates,
    zero_internal_candidates,
)
from order15_mu3.scripts.sat_orthogonal_13 import canonical_second_row, contingency_orbits, contingency_tables


def test_color_14_1_zero_internal_energy_is_arithmetically_impossible() -> None:
    assert zero_internal_candidates()["arithmetic_survivors"] == 0


def test_color_14_1_one_internal_edge_has_ten_arithmetic_candidates() -> None:
    assert one_internal_edge_candidates()["arithmetic_survivors"] == 10


def test_color_14_1_two_internal_edges_are_finitely_enumerated() -> None:
    result = two_internal_edge_candidates()
    assert result["support_types"] == ["two_disjoint_edges", "three_vertex_path"]
    assert result["surviving_determinants"] == sorted(result["surviving_determinants"])


def test_orthogonal_pair_contingency_tables_are_exhaustive() -> None:
    tables = contingency_tables()
    assert len(tables) == 21
    assert all(len(canonical_second_row(table)) == 15 for table in tables)
    assert contingency_orbits() == [
        [0, 5, 20], [1, 4, 6, 10, 18, 19], [2, 3, 11, 14, 15, 17], [7, 9, 16], [8, 12, 13]
    ]
