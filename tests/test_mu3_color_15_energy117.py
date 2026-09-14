from order15_mu3.scripts.verify_color_15_energy117 import (
    MIXED_PARTITIONS,
    NO_ISOLATE_PARTITIONS,
    balanced_unit_support_edge_divisibility,
)
from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions


def test_q117_partition_and_balance_completeness() -> None:
    partitions = set(weight_partitions(13))
    assert len(partitions) == 19
    assert set(MIXED_PARTITIONS) == partitions - {(1,) * 13}
    assert set(NO_ISOLATE_PARTITIONS) == {partition for partition in partitions if len(partition) >= 8}
    assert balanced_unit_support_edge_divisibility()["possible"] is False
