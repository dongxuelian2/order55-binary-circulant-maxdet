from order15_mu3.scripts.verify_q126_k4_four_forest import certificate


def test_q126_k4_four_forest_exact_certificate() -> None:
    result = certificate()
    assert result["graph_type_count"] == 19
    assert result["profile_counts"] == {
        "unit": result["graph_type_count"],
        "weighted": result["graph_type_count"],
        "weighted4": result["graph_type_count"],
    }
    assert result["component_maxima"]["unit"] == 1766431152
    assert result["component_maxima"]["weighted"] == 1628617347
    assert result["component_maxima"]["weighted4"] == 1573793631
    assert result["maximum"] == 256448601928008000
    assert result["weight4_maximum"] == 258141500324775000
    assert result["over_q153"] < 1
