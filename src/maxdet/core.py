"""Core matrix primitives used by the baseline experiments.

The fast routines are intentionally small wrappers around NumPy. They are
appropriate for heuristic ranking, but their floating-point output is never a
mathematical certificate. :func:`exact_determinant` is the certificate path.
"""

from __future__ import annotations

from operator import index as as_index
from typing import Any

import numpy as np
import sympy as sp


def validate_pm1_matrix(matrix: Any) -> np.ndarray:
    """Validate and return a square matrix whose entries are exactly ±1.

    Integer and floating-point arrays are accepted when every value compares
    exactly equal to ``-1`` or ``1``. The returned array is the input view;
    callers that mutate it should make a copy first.
    """

    array = np.asarray(matrix)
    if array.ndim != 2:
        raise ValueError(f"expected a two-dimensional matrix, got ndim={array.ndim}")
    if array.shape[0] != array.shape[1]:
        raise ValueError(f"expected a square matrix, got shape={array.shape}")
    if not (
        np.issubdtype(array.dtype, np.integer)
        or np.issubdtype(array.dtype, np.floating)
    ):
        raise ValueError("matrix entries must use an integer or floating-point dtype")
    if not np.all((array == -1) | (array == 1)):
        raise ValueError("matrix entries must all be exactly -1 or 1")
    return array


def _validated_int8_copy(matrix: Any) -> np.ndarray:
    """Validate a matrix and make a safe mutable ±1 integer copy."""

    return validate_pm1_matrix(matrix).astype(np.int8, copy=True)


def _checked_index(value: Any, size: int, label: str) -> int:
    """Convert an index-like value and provide a consistent bounds error."""

    try:
        position = as_index(value)
    except TypeError as exc:
        raise TypeError(f"{label} index must be an integer") from exc
    if not 0 <= position < size:
        raise IndexError(f"{label} index {position} out of bounds for size {size}")
    return position


def fast_determinant(matrix: Any) -> float:
    """Return a fast floating-point determinant for heuristic scoring only."""

    array = validate_pm1_matrix(matrix).astype(np.float64, copy=False)
    return float(np.linalg.det(array))


def fast_logabs_determinant(matrix: Any) -> float:
    """Return ``log(abs(det(matrix)))`` from NumPy's floating-point ``slogdet``.

    A singular matrix returns ``-inf``. This value is also heuristic and is
    not a mathematical certificate.
    """

    array = validate_pm1_matrix(matrix).astype(np.float64, copy=False)
    _sign, logabsdet = np.linalg.slogdet(array)
    return float(logabsdet)


def exact_determinant(matrix: Any) -> int:
    """Return the determinant as a Python ``int`` using SymPy exact arithmetic."""

    array = validate_pm1_matrix(matrix)
    # Normalizing to int8 makes even an input containing 1.0/-1.0 exact before
    # SymPy sees it, and is safe because validation has already happened.
    integer_matrix = array.astype(np.int8, copy=False).tolist()
    determinant = sp.Matrix(integer_matrix).det(method="domain-ge")
    return int(determinant)


def negate_row(matrix: Any, row_index: int) -> np.ndarray:
    """Return a copy with one row multiplied by -1."""

    result = _validated_int8_copy(matrix)
    row = _checked_index(row_index, result.shape[0], "row")
    result[row, :] = -result[row, :]
    return result


def negate_column(matrix: Any, column_index: int) -> np.ndarray:
    """Return a copy with one column multiplied by -1."""

    result = _validated_int8_copy(matrix)
    column = _checked_index(column_index, result.shape[1], "column")
    result[:, column] = -result[:, column]
    return result


def swap_rows(matrix: Any, first: int, second: int) -> np.ndarray:
    """Return a copy with two rows exchanged."""

    result = _validated_int8_copy(matrix)
    first_row = _checked_index(first, result.shape[0], "row")
    second_row = _checked_index(second, result.shape[0], "row")
    result[[first_row, second_row], :] = result[[second_row, first_row], :]
    return result


def swap_columns(matrix: Any, first: int, second: int) -> np.ndarray:
    """Return a copy with two columns exchanged."""

    result = _validated_int8_copy(matrix)
    first_column = _checked_index(first, result.shape[1], "column")
    second_column = _checked_index(second, result.shape[1], "column")
    result[:, [first_column, second_column]] = result[:, [second_column, first_column]]
    return result


def normalize_pm1_matrix(matrix: Any) -> np.ndarray:
    """Normalize signs so the first row and first column are all +1.

    This uses only row and column sign flips, so the absolute determinant is
    unchanged. No canonical labeling or equivalence classification is done.
    """

    result = _validated_int8_copy(matrix)
    result *= result[:, [0]]
    result *= result[[0], :]
    return result


def sylvester_hadamard(order: int) -> np.ndarray:
    """Construct a Sylvester Hadamard matrix of order ``1, 2, 4, ...``."""

    try:
        size = as_index(order)
    except TypeError as exc:
        raise TypeError("Hadamard order must be an integer") from exc
    if size < 1 or size & (size - 1):
        raise ValueError("Sylvester Hadamard order must be a positive power of two")

    matrix = np.array([[1]], dtype=np.int8)
    while matrix.shape[0] < size:
        matrix = np.concatenate(
            [
                np.concatenate([matrix, matrix], axis=1),
                np.concatenate([matrix, -matrix], axis=1),
            ],
            axis=0,
        )
    return matrix


def random_pm1_matrix(
    n: int,
    rng: np.random.Generator | int | None,
) -> np.ndarray:
    """Generate an ``n × n`` random ±1 matrix from an explicit RNG or seed."""

    try:
        size = as_index(n)
    except TypeError as exc:
        raise TypeError("matrix order n must be an integer") from exc
    if size < 1:
        raise ValueError("matrix order n must be positive")

    generator = rng if isinstance(rng, np.random.Generator) else np.random.default_rng(rng)
    bits = generator.integers(0, 2, size=(size, size), dtype=np.int8)
    return (2 * bits - 1).astype(np.int8, copy=False)
