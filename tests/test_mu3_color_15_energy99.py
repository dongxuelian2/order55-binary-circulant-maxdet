from order15_mu3.scripts.verify_color_15_energy99 import Q99_ALL_SAME_UPPER, certificate


def test_q99_exact_all_same_boundary_certificate() -> None:
    result = certificate()
    assert result["isolated_rank_pair_survivors"] == []
    assert result["all_same_color_upper"] == Q99_ALL_SAME_UPPER
