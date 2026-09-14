from order15_mu3.scripts.verify_aggregate_upper import certificate


NEW_UPPER = 287677100533178196
OLD_UPPER = 287710229992756239


def test_exact_branch_aggregate_after_tail_closure() -> None:
    result = certificate()
    assert result["closed_post_q150_shells"] == [153, 159, 162, 168]
    assert result["counterexample_energy_bound"] == "Q <= 150"
    assert result["maximizing_source"] == "Q=150 exact boundary"
    assert result["rigorous_integer_upper"] == NEW_UPPER
    assert result["non_norm_multiples_skipped"] == 6
    assert result["rigorous_integer_upper"] < OLD_UPPER
    assert result["normalized_lower_percent"] < result["normalized_upper_percent"]
    assert result["normalized_upper_percent"] < 81.053
