from order15_mu3.scripts.verify_color_14_1_twin_bound import certificates


def test_q60_and_q69_twin_certificates_are_negative() -> None:
    rows = certificates()
    assert {row["total_energy"] for row in rows} == {60, 69}
    assert all(row["forced_gram_eigenvalue_upper"] < 0 for row in rows)
