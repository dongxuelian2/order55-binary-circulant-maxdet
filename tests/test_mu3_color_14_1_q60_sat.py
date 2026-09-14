import pytest

from order15_mu3.scripts.sat_color_14_1_q60 import build, gram_value


def test_q60_simultaneous_gram_model_builds() -> None:
    cnf, pool, _x = build("minus", "minus")
    assert cnf.clauses and pool.top > 0
    assert gram_value(0, 14, "minus").norm() == 3
    assert gram_value(0, 1, "minus").norm() == 9


def test_q60_rejects_different_determinant_classes() -> None:
    with pytest.raises(ValueError):
        build("minus", "plus")
