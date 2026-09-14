from order15_mu3.scripts.verify_color_15_energy108 import MIXED_PARTITIONS
from order15_mu3.scripts.verify_color_15_energy90 import weight_partitions


def test_q108_weight_partition_completeness() -> None:
    partitions = set(weight_partitions(12))
    assert len(partitions) == 16
    assert set(MIXED_PARTITIONS) == partitions - {(1,) * 12}
