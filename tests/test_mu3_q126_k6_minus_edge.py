from order15_mu3.scripts.verify_q126_k6_minus_edge import certificate


def test_q126_k6_minus_edge_exact_phase_certificate() -> None:
    result = certificate()
    assert result["switching_normalized_phase_assignments"] == 6**9
    assert result["over_q153"] < 1
