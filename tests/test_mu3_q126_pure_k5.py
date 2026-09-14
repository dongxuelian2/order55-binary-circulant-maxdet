from order15_mu3.scripts.verify_q126_pure_k5 import certificate


def test_q126_pure_k5_exact_certificate() -> None:
    result = certificate()
    assert result["phase_profile_count"] == 1
    assert result["component_maximum"] == 559872
    assert result["maxima"]["joined0"] == 274211883004723200
    assert result["leaf_component_maximum"] == 8024832
    assert result["maxima"]["joined1_leaf"] == 272942383546368000
    assert result["over_q153"] < 1
