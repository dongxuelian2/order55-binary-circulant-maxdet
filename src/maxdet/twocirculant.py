"""Bordered two-circulant constructions used in recent maxdet records."""

from __future__ import annotations

from typing import Iterable

import numpy as np


def parse_sign_word(word: str) -> np.ndarray:
    """Parse a compact ``+``/``-`` word as an int8 vector."""
    if not word or any(character not in "+-" for character in word):
        raise ValueError("sign word must be a nonempty string containing only '+' and '-'")
    return np.fromiter((1 if character == "+" else -1 for character in word), dtype=np.int8)


def format_sign_word(values: Iterable[int]) -> str:
    """Format a one-dimensional ±1 iterable as a compact sign word."""
    array = np.asarray(list(values))
    if array.ndim != 1 or array.size == 0 or not np.all((array == -1) | (array == 1)):
        raise ValueError("values must be a nonempty one-dimensional ±1 sequence")
    return "".join("+" if value == 1 else "-" for value in array)


def circulant_from_first_row(row: Iterable[int]) -> np.ndarray:
    """Return the circulant convention used by the Caltech record verifier."""
    values = np.asarray(list(row), dtype=np.int8)
    if values.ndim != 1 or values.size == 0 or not np.all((values == -1) | (values == 1)):
        raise ValueError("row must be a nonempty one-dimensional ±1 sequence")
    indices = np.arange(values.size)
    return values[(indices[None, :] - indices[:, None]) % values.size]


def bordered_two_circulant(a: Iterable[int], b: Iterable[int], matrix_type: int) -> np.ndarray:
    """Construct a ``(2m+1) × (2m+1)`` ±1 matrix of public type 0, 1, 2, or 3."""
    first = np.asarray(list(a), dtype=np.int8)
    second = np.asarray(list(b), dtype=np.int8)
    if first.ndim != 1 or second.ndim != 1 or first.size == 0 or first.shape != second.shape:
        raise ValueError("a and b must be nonempty one-dimensional words of equal length")
    if not np.all((first == -1) | (first == 1)) or not np.all((second == -1) | (second == 1)):
        raise ValueError("a and b entries must all be ±1")
    if matrix_type not in (0, 1, 2, 3):
        raise ValueError("matrix_type must be one of 0, 1, 2, 3")
    A = circulant_from_first_row(first)
    B = circulant_from_first_row(second)
    m = first.size
    column = np.ones((m, 1), dtype=np.int8)
    scalar = np.ones((1, 1), dtype=np.int8)
    if matrix_type == 0:
        return np.block([[-scalar, column.T, -column.T], [column, A, B], [-column, B.T, -A.T]])
    if matrix_type == 1:
        return np.block([[scalar, column.T, -column.T], [column, A, B], [-column, B.T, -A.T]])
    if matrix_type == 2:
        return np.block([[A, B, -column], [B, A.T, column], [column.T, -column.T, scalar]])
    return np.block([[A, B, -column], [B.T, -A.T, column], [-column.T, -column.T, -scalar]])


def spectral_logabsdet(a: Iterable[int], b: Iterable[int], matrix_type: int) -> float:
    """Heuristically score type-1 or type-3 layouts by Fourier block diagonalization."""
    first = np.asarray(list(a), dtype=np.float64)
    second = np.asarray(list(b), dtype=np.float64)
    if first.ndim != 1 or first.shape != second.shape or first.size % 2 != 1:
        raise ValueError("a and b must be equal-length odd one-dimensional words")
    if matrix_type not in (1, 3):
        raise ValueError("spectral scorer currently supports matrix types 1 and 3")
    if not np.all((first == -1) | (first == 1)) or not np.all((second == -1) | (second == 1)):
        raise ValueError("a and b entries must all be ±1")
    m = first.size
    alpha = float(first.sum())
    beta = float(second.sum())
    if matrix_type == 1:
        zero_block = np.array([[1.0, m, -m], [1.0, alpha, beta], [-1.0, beta, -alpha]])
    else:
        zero_block = np.array([[alpha, beta, -1.0], [beta, -alpha, 1.0], [-m, -m, -1.0]])
    sign, result = np.linalg.slogdet(zero_block)
    if sign == 0:
        return float("-inf")
    spectrum_a = np.fft.fft(first)
    spectrum_b = np.fft.fft(second)
    for frequency in range(1, (m + 1) // 2):
        magnitude = abs(spectrum_a[frequency]) ** 2 + abs(spectrum_b[frequency]) ** 2
        if magnitude == 0:
            return float("-inf")
        result += 2.0 * np.log(magnitude)
    return float(result)
