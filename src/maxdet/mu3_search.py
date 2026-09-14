"""Floating-point discovery helpers for dephased third-root matrices.

These routines rank candidates only. Any record must be certified by the exact
routines in :mod:`maxdet.mu3`.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np


COMPLEX_ROOTS = np.exp(2j * np.pi * np.arange(3) / 3)


def complex_matrix(exponents: np.ndarray) -> np.ndarray:
    return COMPLEX_ROOTS[np.asarray(exponents, dtype=np.int8) % 3]


def random_dephased(order: int, rng: np.random.Generator) -> np.ndarray:
    x = rng.integers(0, 3, size=(order, order), dtype=np.int8)
    x[0, :] = 0
    x[:, 0] = 0
    return x


@dataclass
class LocalResult:
    exponents: np.ndarray
    logabsdet: float
    accepted_moves: int


def single_entry_ascent(exponents: np.ndarray, tolerance: float = 1e-12) -> LocalResult:
    """Reach a strict single-entry local maximum using inverse updates."""

    x = np.asarray(exponents, dtype=np.int8).copy()
    h = complex_matrix(x)
    sign, logabs = np.linalg.slogdet(h)
    if sign == 0:
        raise np.linalg.LinAlgError("singular starting matrix")
    inverse = np.linalg.inv(h)
    accepted = 0
    while True:
        best: tuple[float, int, int, int, complex] | None = None
        for i in range(1, x.shape[0]):
            for j in range(1, x.shape[1]):
                old = int(x[i, j])
                for new in ((old + 1) % 3, (old + 2) % 3):
                    delta = COMPLEX_ROOTS[new] - COMPLEX_ROOTS[old]
                    ratio = 1.0 + delta * inverse[j, i]
                    ratio_abs = abs(ratio)
                    gain = float(np.log(ratio_abs)) if ratio_abs else float("-inf")
                    if gain > tolerance and (best is None or gain > best[0]):
                        best = (gain, i, j, new, ratio)
        if best is None:
            return LocalResult(x, float(logabs), accepted)
        gain, i, j, new, ratio = best
        old_value = h[i, j]
        new_value = COMPLEX_ROOTS[new]
        delta = new_value - old_value
        inverse -= (delta / ratio) * np.outer(inverse[:, i], inverse[j, :])
        h[i, j] = new_value
        x[i, j] = new
        logabs += gain
        accepted += 1
        if accepted % 100 == 0:
            sign, logabs = np.linalg.slogdet(h)
            if sign == 0:
                raise np.linalg.LinAlgError("numerical singularity during ascent")
            inverse = np.linalg.inv(h)


def perturb(exponents: np.ndarray, changes: int, rng: np.random.Generator) -> np.ndarray:
    result = np.asarray(exponents, dtype=np.int8).copy()
    order = result.shape[0]
    positions = rng.choice((order - 1) ** 2, size=min(changes, (order - 1) ** 2), replace=False)
    for position in positions:
        i, j = divmod(int(position), order - 1)
        i += 1
        j += 1
        result[i, j] = (int(result[i, j]) + int(rng.integers(1, 3))) % 3
    return result


def best_dephased_row(coefficients: np.ndarray, tolerance: float = 1e-11) -> tuple[np.ndarray, complex]:
    """Maximize ``abs(sum(root[x[j]]*coefficients[j]))`` with ``x[0]=0``.

    The objective is the support function of a Minkowski sum of triangles.
    Between angular tie events every coordinate has a fixed maximizing phase;
    at a simultaneous event we explicitly enumerate the tied choices.
    """

    coefficients = np.asarray(coefficients, dtype=np.complex128)
    if coefficients.ndim != 1 or coefficients.size == 0:
        raise ValueError("coefficients must be a nonempty vector")
    boundaries: list[float] = [0.0]
    for coefficient in coefficients[1:]:
        if abs(coefficient) <= tolerance:
            continue
        for first in range(3):
            for second in range(first + 1, 3):
                angle = np.angle((COMPLEX_ROOTS[first] - COMPLEX_ROOTS[second]) * coefficient)
                boundaries.append(float((angle + np.pi / 2) % (2 * np.pi)))
                boundaries.append(float((angle - np.pi / 2) % (2 * np.pi)))
    boundaries = sorted(boundaries)
    unique: list[float] = []
    for angle in boundaries:
        if not unique or abs(angle - unique[-1]) > tolerance:
            unique.append(angle)

    candidate_angles = unique[:]
    for index, left in enumerate(unique):
        right = unique[(index + 1) % len(unique)]
        if index + 1 == len(unique):
            right += 2 * np.pi
        candidate_angles.append(float(((left + right) / 2) % (2 * np.pi)))

    best_exponents = np.zeros(coefficients.size, dtype=np.int8)
    best_sum = complex(np.sum(coefficients))
    best_abs = abs(best_sum)
    for theta in candidate_angles:
        direction = np.exp(-1j * theta)
        choices: list[list[int]] = [[0]]
        combinations = 1
        for coefficient in coefficients[1:]:
            projections = np.real(direction * COMPLEX_ROOTS * coefficient)
            maximum = float(np.max(projections))
            tied = [index for index, value in enumerate(projections) if maximum - float(value) <= tolerance]
            if abs(coefficient) <= tolerance:
                tied = [0]
            choices.append(tied)
            combinations *= len(tied)
        # At most 3^14 in a fully degenerate case; zero coefficients were fixed
        # above, and genuine large simultaneous ties are rare at order 15.
        if combinations > 200_000:
            choices = [[options[0]] for options in choices]
        for selection in product(*choices):
            exponents = np.asarray(selection, dtype=np.int8)
            value = complex(np.dot(COMPLEX_ROOTS[exponents], coefficients))
            magnitude = abs(value)
            if magnitude > best_abs + tolerance:
                best_abs = magnitude
                best_sum = value
                best_exponents = exponents
    return best_exponents, best_sum


def alternating_row_column_ascent(
    exponents: np.ndarray,
    max_passes: int = 100,
    tolerance: float = 1e-10,
) -> LocalResult:
    """Alternate globally optimal whole-row and whole-column replacements."""

    x = np.asarray(exponents, dtype=np.int8).copy()
    accepted = 0
    for _ in range(max_passes):
        changed = False
        for transpose in (False, True):
            view = x.T if transpose else x
            for row in range(1, view.shape[0]):
                h = complex_matrix(view)
                sign, logabs = np.linalg.slogdet(h)
                if sign == 0:
                    raise np.linalg.LinAlgError("singular matrix during row/column ascent")
                coefficients = np.linalg.inv(h)[:, row]
                proposal, ratio = best_dephased_row(coefficients)
                if abs(ratio) > 1 + tolerance and not np.array_equal(proposal, view[row]):
                    view[row, :] = proposal
                    accepted += 1
                    changed = True
        if not changed:
            break
    _sign, logabs = np.linalg.slogdet(complex_matrix(x))
    return LocalResult(x, float(logabs), accepted)
