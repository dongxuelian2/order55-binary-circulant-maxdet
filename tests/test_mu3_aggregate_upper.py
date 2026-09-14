from order15_mu3.scripts.verify_aggregate_upper import certificate


def test_exact_branch_aggregate_improves_global_upper() -> None:
    result = certificate()
    assert result["rigorous_integer_upper"] < 293695860429200709
    assert result["squared_gap_factor"] < 1.057
