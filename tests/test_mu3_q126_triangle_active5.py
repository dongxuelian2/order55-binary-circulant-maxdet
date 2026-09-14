from order15_mu3.scripts.verify_q126_triangle_active5 import certificate


def test_q126_triangle_active5_exact_certificate() -> None:
    result = certificate()
    assert len(result["maxima"]) == 3
    assert result["graph_type_counts"] == {
        "edges10_component5": 3,
        "edges11_component6": 4,
        "edges12_component7": 3,
    }
    assert result["maximum"] == 255359816048148480
    assert result["over_q153"] < 1
