import numpy as np

from itertools import product

from maxdet.mu3_search import (
    COMPLEX_ROOTS,
    alternating_row_column_ascent,
    best_dephased_row,
    complex_matrix,
    random_dephased,
    single_entry_ascent,
)


def test_rank_one_scores_and_local_ascent() -> None:
    rng = np.random.default_rng(7)
    while True:
        x = random_dephased(5, rng)
        h = complex_matrix(x)
        if abs(np.linalg.det(h)) > 1e-8:
            break
    before = abs(np.linalg.det(h))
    result = single_entry_ascent(x)
    after = abs(np.linalg.det(complex_matrix(result.exponents)))
    assert after + 1e-8 >= before
    inverse = np.linalg.inv(complex_matrix(result.exponents))
    for i in range(1, 5):
        for j in range(1, 5):
            old = int(result.exponents[i, j])
            for new in ((old + 1) % 3, (old + 2) % 3):
                delta = COMPLEX_ROOTS[new] - COMPLEX_ROOTS[old]
                assert abs(1 + delta * inverse[j, i]) <= 1 + 1e-9


def test_angular_row_optimizer_matches_brute_force() -> None:
    rng = np.random.default_rng(19)
    for _ in range(20):
        coefficients = rng.normal(size=6) + 1j * rng.normal(size=6)
        exponents, value = best_dephased_row(coefficients)
        brute = max(
            abs(np.dot(COMPLEX_ROOTS[np.asarray((0, *tail), dtype=np.int8)], coefficients))
            for tail in product(range(3), repeat=5)
        )
        assert exponents[0] == 0
        assert abs(abs(value) - brute) < 1e-9


def test_alternating_row_column_ascent_is_monotone() -> None:
    rng = np.random.default_rng(23)
    while True:
        x = random_dephased(6, rng)
        before = abs(np.linalg.det(complex_matrix(x)))
        if before > 1e-8:
            break
    result = alternating_row_column_ascent(x)
    after = abs(np.linalg.det(complex_matrix(result.exponents)))
    assert after + 1e-8 >= before
    assert np.all(result.exponents[0, :] == 0)
    assert np.all(result.exponents[:, 0] == 0)
