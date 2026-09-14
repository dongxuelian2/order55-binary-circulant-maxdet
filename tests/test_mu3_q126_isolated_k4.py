from order15_mu3.scripts.verify_q126_isolated_k4 import certificate


def test_q126_isolated_k4_exact_certificate() -> None:
    result = certificate()
    assert set(result["component_maxima"]) == {"unit", "weighted3", "weighted4"}
    assert set(result["maxima"]) == {"edges11_weight4", "edges12_weight3"}
    assert result["component_maxima"]["unit"] == 41472
    assert result["component_maxima"]["weighted3"] == 38232
    assert result["component_maxima"]["weighted4"] == 36936
    assert result["maxima"] == {
        "edges11_weight4": 260501288854487040,
        "edges12_weight3": 258406614748200960,
    }
    assert result["over_q153"] < 1
