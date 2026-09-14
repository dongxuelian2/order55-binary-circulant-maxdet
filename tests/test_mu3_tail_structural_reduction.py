from order15_mu3.scripts.verify_tail_structural_reduction import certificate


def test_tail_structural_reduction_closes_q162_and_q168() -> None:
    result = certificate()

    assert result["excluded_tail_shells"] == [162, 168]
    assert result["remaining_tail_shells"] == [153, 159]

    assert result["q162"]["all_same_over_record"] < 1
    assert all(row["excluded"] for row in result["q162"]["size13_rows"])
    assert all(row["excluded"] for row in result["q168"]["size13_rows"])
    assert all(row["excluded"] for row in result["q168"]["size14_rows"])


def test_tail_residual_frontier_is_explicit() -> None:
    result = certificate()

    q153 = result["q153"]
    assert q153["residual_branch_count"] == 5
    assert {
        (row["internal"], row["cross"], row["outside_norm"])
        for row in q153["residual_size13_allocations"]
    } == {
        (27, 78, 48),
        (27, 87, 39),
        (36, 78, 39),
        (54, 78, 21),
    }

    q159 = result["q159"]
    assert q159["residual_branch_count"] == 4
    assert {
        (row["internal"], row["cross"], row["outside_norm"])
        for row in q159["residual_size13_allocations"]
    } == {(54, 78, 27)}
    assert {row["internal"] for row in q159["residual_size14_allocations"]} == {
        99,
        108,
        117,
    }
