from order15_mu3.scripts.verify_support_independence import independence_number


def test_low_energy_support_independence_numbers() -> None:
    assert independence_number(14, ((0, 1), (2, 3))) == 12
    assert independence_number(14, ((0, 1), (1, 2))) == 13
    assert independence_number(14, ((0, 1), (0, 2), (0, 3))) == 13
