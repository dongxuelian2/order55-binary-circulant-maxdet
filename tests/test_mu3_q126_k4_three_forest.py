from order15_mu3.scripts.verify_q126_k4_three_forest import certificate


def test_q126_k4_three_forest_exact_certificate() -> None:
    result = certificate()
    assert result["graph_type_count"] == result["profile_counts"]["unit"]
    assert result["graph_type_count"] == result["profile_counts"]["weighted"]
    assert result["graph_type_count"] == 7
    assert result["component_maxima"] == {"unit": 122891904, "weighted": 113304096}
    assert result["maximum"] == 256914952884633600
    assert result["over_q153"] < 1
