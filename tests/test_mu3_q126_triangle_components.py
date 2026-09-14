from order15_mu3.scripts.verify_q126_triangle_components import certificate


def test_q126_triangle_component_exact_certificate() -> None:
    result = certificate()
    assert result["graph_type_counts"] == {
        "active10_unicyclic": 299,
        "active9_unicyclic": 117,
        "active8_bicyclic": 173,
        "active7_bicyclic": 51,
    }
    assert result["maxima"] == {
        "active10_unicyclic": 246462947301769920,
        "active9_unicyclic": 246911098793504256,
        "active8_bicyclic": 251108976455562240,
        "active7_bicyclic": 251565599541989376,
    }
    assert result["over_q153"] < 1
