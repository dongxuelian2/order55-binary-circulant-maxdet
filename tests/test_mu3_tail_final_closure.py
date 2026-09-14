from order15_mu3.scripts.verify_tail_final_closure import certificate


def test_all_genuine_post_q150_shells_are_closed() -> None:
    result = certificate()
    assert result["closed_tail_shells"] == [153, 159, 162, 168]
    assert result["new_counterexample_energy_bound"] == "Q <= 150"
    assert all(value < 1 for value in result["q153_size13_improved_over_record"].values())
    assert all(value < 1 for value in result["q159_size13_improved_over_record"].values())


def test_q153_all_same_boundary_is_below_record() -> None:
    row = certificate()["q153_all_same"]
    assert row["clique_at_most_five_over_record"] < 1
    assert row["eight_or_more_isolates_over_record"] < 1
    assert row["five_isolates_over_record"] < 1
    assert row["six_isolates_over_record"] < 1
    assert row["seven_isolates_adjacent_over_record"] < 1
    assert row["seven_isolates_cross_over_record"] < 1


def test_q159_size14_component_reduction_closes_boundary() -> None:
    row = certificate()["q159_size14"]
    assert {entry["internal"] for entry in row["rows"]} == {99, 108, 117}
    assert row["rows"][0]["generic_residual_configurations"] == []
    assert row["rows"][1]["generic_residual_configurations"] == []
    assert row["rows"][2]["generic_residual_configurations"] == [[[14, 13]]]
    assert row["connected_e117_tree_boundary"]["unlabeled_trees"] == 3159
    assert row["connected_e117_tree_boundary"]["maximizer"] == "P14"
    assert row["connected_e117_tree_boundary"]["over_record"] < 1
