from order15_mu3.scripts.verify_q87_q90_second_obstruction import NEXT, TARGET, certificate


def test_second_energy9_determinant_is_excluded() -> None:
    result = certificate()
    assert result["excluded_determinant"] == TARGET
    assert result["next_arithmetic_energy9_upper"] == NEXT
    assert all(
        row["result"] == "unsat"
        for case in result["color_cases"].values()
        for row in case["pair_results"]
    )
