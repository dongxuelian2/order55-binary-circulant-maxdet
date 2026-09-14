from order15_mu3.scripts.verify_q141_energy36_minimal_cross import certificate


def test_q141_energy36_minimal_cross_exact_certificate() -> None:
    result = certificate()
    assert result["profiles"]["K2+K3"]["maximum"] == 248427409860000000
    assert result["over_q153"] < 1
