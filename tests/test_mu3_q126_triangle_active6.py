from order15_mu3.scripts.verify_q126_triangle_active6 import certificate


def test_q126_triangle_active6_exact_certificate() -> None:
    result = certificate()
    assert set(result["maxima"]) == {
        "edges10_component6",
        "edges11_component7",
        "edges12_component8",
    }
    assert result["graph_type_counts"] == {
        "edges10_component6": 7,
        "edges11_component7": 15,
        "edges12_component8": 17,
    }
    assert result["maximum"] == 255169391129395200
    assert result["over_q153"] < 1
