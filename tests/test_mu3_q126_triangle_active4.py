from order15_mu3.scripts.verify_q126_triangle_active4 import certificate


def test_q126_triangle_active4_exact_certificate() -> None:
    result = certificate()
    assert set(result["maxima"]) == {"edges9_component4", "edges10_component5"}
    assert result["graph_type_counts"] == {"edges9_component4": 1, "edges10_component5": 1}
    assert result["maximum"] == 254788541291888640
    assert result["over_q153"] < 1
