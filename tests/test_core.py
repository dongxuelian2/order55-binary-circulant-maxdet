from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys

import numpy as np
import pytest

from maxdet import (
    exact_determinant,
    fast_determinant,
    fast_logabs_determinant,
    negate_column,
    negate_row,
    normalize_pm1_matrix,
    random_pm1_matrix,
    run_random_baseline,
    swap_columns,
    swap_rows,
    sylvester_hadamard,
    validate_pm1_matrix,
)


def _sample_matrix() -> np.ndarray:
    return np.array(
        [
            [1, -1, 1, 1],
            [-1, 1, 1, -1],
            [1, 1, -1, 1],
            [1, -1, -1, -1],
        ],
        dtype=np.int8,
    )


def test_pm1_validation_accepts_square_pm1_matrix() -> None:
    matrix = validate_pm1_matrix([[1, -1], [-1, 1]])
    assert matrix.shape == (2, 2)


@pytest.mark.parametrize(
    "matrix",
    [
        [1, -1, 1],
        [[1, 0], [-1, 1]],
        [[1, -1, 1], [-1, 1, -1]],
        [[True, False], [False, True]],
    ],
)
def test_pm1_validation_rejects_invalid_matrix(matrix) -> None:
    with pytest.raises(ValueError):
        validate_pm1_matrix(matrix)


def test_exact_determinant_hand_check() -> None:
    matrix = np.array([[1, 1], [1, -1]], dtype=np.int8)
    assert exact_determinant(matrix) == -2


def test_fast_and_exact_determinants_agree_on_small_matrices() -> None:
    rng = np.random.default_rng(17)
    for n in range(1, 6):
        for _ in range(4):
            matrix = random_pm1_matrix(n, rng)
            exact = exact_determinant(matrix)
            assert fast_determinant(matrix) == pytest.approx(exact)
            if exact == 0:
                assert fast_logabs_determinant(matrix) == -np.inf
            else:
                assert fast_logabs_determinant(matrix) == pytest.approx(np.log(abs(exact)))


@pytest.mark.parametrize(
    "transform",
    [
        lambda matrix: negate_row(matrix, 1),
        lambda matrix: negate_column(matrix, 2),
        lambda matrix: swap_rows(matrix, 0, 1),
        lambda matrix: swap_columns(matrix, 0, 2),
    ],
)
def test_basic_symmetries_preserve_absolute_determinant(transform) -> None:
    matrix = _sample_matrix()
    transformed = transform(matrix)
    assert np.array_equal(matrix, _sample_matrix())
    assert abs(exact_determinant(transformed)) == abs(exact_determinant(matrix))


def test_normalization_sets_first_row_and_column_and_preserves_absolute_det() -> None:
    matrix = _sample_matrix()
    normalized = normalize_pm1_matrix(matrix)
    assert np.all(normalized[0, :] == 1)
    assert np.all(normalized[:, 0] == 1)
    assert abs(exact_determinant(normalized)) == abs(exact_determinant(matrix))


@pytest.mark.parametrize("order", [1, 2, 4, 8])
def test_sylvester_hadamard_fixture(order: int) -> None:
    matrix = sylvester_hadamard(order)
    expected = 1 if order == 1 else order ** (order // 2)
    assert matrix.shape == (order, order)
    assert np.all((matrix == -1) | (matrix == 1))
    assert abs(exact_determinant(matrix)) == expected


def _normalized_candidates(order: int):
    inner = order - 1
    free_entries = inner * inner
    for mask in range(1 << free_entries):
        matrix = np.ones((order, order), dtype=np.int8)
        for position in range(free_entries):
            row, column = divmod(position, inner) if inner else (0, 0)
            if (mask >> position) & 1:
                matrix[row + 1, column + 1] = -1
        yield matrix


@pytest.mark.parametrize("order,expected", [(1, 1), (2, 2), (3, 4), (4, 16)])
def test_known_maxima_by_tiny_normalized_enumeration(order: int, expected: int) -> None:
    best_score = -1.0
    best_matrix = None
    for matrix in _normalized_candidates(order):
        score = abs(fast_determinant(matrix))
        if score > best_score:
            best_score = score
            best_matrix = matrix
    assert best_matrix is not None
    assert abs(exact_determinant(best_matrix)) == expected


def test_rng_reproducibility() -> None:
    first = random_pm1_matrix(6, np.random.default_rng(1234))
    second = random_pm1_matrix(6, np.random.default_rng(1234))
    assert np.array_equal(first, second)

    first_run = run_random_baseline(n=5, samples=12, seed=9)
    second_run = run_random_baseline(n=5, samples=12, seed=9)
    assert first_run.best_fast_score == second_run.best_fast_score
    assert first_run.exact_determinant == second_run.exact_determinant
    assert np.array_equal(first_run.matrix, second_run.matrix)


def test_random_baseline_cli_smoke() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    script = repo_root / "scripts" / "random_baseline.py"
    environment = os.environ.copy()
    source_path = str(repo_root / "src")
    environment["PYTHONPATH"] = source_path + os.pathsep + environment.get("PYTHONPATH", "")
    completed = subprocess.run(
        [sys.executable, str(script), "--n", "4", "--samples", "8", "--seed", "0"],
        cwd=repo_root,
        check=True,
        capture_output=True,
        text=True,
        env=environment,
    )
    assert "n: 4" in completed.stdout
    assert "samples: 8" in completed.stdout
    assert "seed: 0" in completed.stdout
    assert "best fast score:" in completed.stdout
    assert "exact determinant:" in completed.stdout
    assert "matrix:" in completed.stdout
