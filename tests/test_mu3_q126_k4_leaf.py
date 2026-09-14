from order15_mu3.scripts.verify_q126_k4_leaf import certificate


def test_q126_k4_leaf_exact_phase_certificate() -> None:
    result = certificate()
    assert result["maximum"] == 257873424975691776
    assert result["over_q153"] < 1
