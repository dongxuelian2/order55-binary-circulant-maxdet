import numpy as np

from maxdet.mu3 import determinant, exponent_matrix
from order15_mu3.scripts.circulant_scan import circulant, digits_base_three


def test_base_three_enumeration_and_circulant_exact_value() -> None:
    rows = digits_base_three(0, 10)
    assert len({tuple(row) for row in rows}) == 10
    assert np.all(rows[:, 0] == 0)
    first = np.asarray([0, 0, 1], dtype=np.int8)
    value = determinant(exponent_matrix(circulant(first)))
    assert value.norm() == 27
