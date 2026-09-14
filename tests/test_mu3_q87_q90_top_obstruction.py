from order15_mu3.scripts.verify_q87_q90_top_obstruction import NEXT, TOP, certificate


def test_top_energy9_determinant_is_excluded() -> None:
    result = certificate()
    assert result["excluded_determinant"] == TOP
    assert result["next_arithmetic_energy9_upper"] == NEXT
    assert all(
        all(row["result"] == "unsat" for row in case["pair_results"])
        for case in result["color_cases"].values()
    )
