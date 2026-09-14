from order15_mu3.scripts.verify_color_15_energy90 import (
    SIZE13_Q90_MAXIMUM,
    certificate,
    weight_partitions,
)


def test_q90_weight_partition_completeness() -> None:
    partitions = set(weight_partitions(10))
    assert len(partitions) == 11
    assert all(sum(partition) == 10 for partition in partitions)


def test_q90_exact_boundary_certificate() -> None:
    result = certificate()
    assert result["isolated_rank_pair_survivors"] == []
    assert result["all_Q90_upper"] == SIZE13_Q90_MAXIMUM
