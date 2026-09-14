from order15_mu3.scripts.verify_size13_exact_loss import EXPECTED_LOCAL_LOSSES, certificate


def test_size13_exact_local_losses() -> None:
    result = certificate()
    assert set(result["local_minimum_losses"]) == set(EXPECTED_LOCAL_LOSSES)
    assert max(result["energy18_bounds_over_q135"].values()) < 1
    assert max(result["energy27_bounds_over_q135"].values()) < 1
