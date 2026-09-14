from order15_mu3.scripts.verify_q126_pure_k4_six_forest import certificate


def test_q126_pure_k4_six_forest_exact_certificate() -> None:
    result = certificate()
    assert result["graph_type_count"] == result["profile_count"]
    assert result["graph_type_count"] == 124
    assert result["component_maximum"] == 364958721792
    assert result["maximum"] == 255412711858913280
    assert result["over_q153"] < 1
