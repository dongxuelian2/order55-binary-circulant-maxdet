"""Exact phase enumeration for the dense Q=126 K6-minus-edge support."""

from __future__ import annotations

import json
from functools import lru_cache
from itertools import permutations

import numpy as np

try:
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_trace_stability import stationary_upper


ORDER = 6
MISSING_EDGE = (4, 5)
TREE_EDGES = tuple((0, vertex) for vertex in range(1, ORDER))
CYCLE_EDGES = tuple(
    (i, j)
    for i in range(ORDER)
    for j in range(i + 1, ORDER)
    if (i, j) not in TREE_EDGES and (i, j) != MISSING_EDGE
)


def _permutation_sign(values: tuple[int, ...]) -> int:
    inversions = sum(
        values[i] > values[j]
        for i in range(len(values))
        for j in range(i + 1, len(values))
    )
    return -1 if inversions % 2 else 1


def determinant_polynomial() -> dict[tuple[int, ...], int]:
    """Laurent polynomial for det(5I+A), with five tree phases gauged to 1."""

    edge_index = {edge: index for index, edge in enumerate(CYCLE_EDGES)}
    coefficients: dict[tuple[int, ...], int] = {}
    for permutation in permutations(range(ORDER)):
        coefficient = _permutation_sign(permutation)
        exponents = [0] * len(CYCLE_EDGES)
        for row, column in enumerate(permutation):
            if row == column:
                coefficient *= 5
                continue
            edge = tuple(sorted((row, column)))
            if edge == MISSING_EDGE:
                coefficient = 0
                break
            if edge in edge_index:
                exponents[edge_index[edge]] += 1 if row < column else -1
        if coefficient:
            key = tuple(exponent % 6 for exponent in exponents)
            coefficients[key] = coefficients.get(key, 0) + coefficient
    coefficients = {key: value for key, value in coefficients.items() if value}
    assert len(CYCLE_EDGES) == 9
    assert len(coefficients) == 291
    return coefficients


# Powers of the primitive sixth root -omega^2=1+omega, stored as a+b*omega.
ROOTS = ((1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1))


def _transform_axis(a: np.ndarray, b: np.ndarray, axis: int) -> tuple[np.ndarray, np.ndarray]:
    """Exact length-six Fourier transform over Eisenstein integer pairs."""

    source_a = np.moveaxis(a, axis, 0).reshape(6, -1)
    source_b = np.moveaxis(b, axis, 0).reshape(6, -1)
    target_a = np.zeros_like(source_a)
    target_b = np.zeros_like(source_b)
    for frequency in range(6):
        for exponent in range(6):
            root_a, root_b = ROOTS[(frequency * exponent) % 6]
            left = source_a[exponent]
            right = source_b[exponent]
            target_a[frequency] += left * root_a - right * root_b
            target_b[frequency] += left * root_b + right * root_a - right * root_b
    shape = np.moveaxis(a, axis, 0).shape
    return (
        np.moveaxis(target_a.reshape(shape), 0, axis),
        np.moveaxis(target_b.reshape(shape), 0, axis),
    )


@lru_cache(maxsize=1)
def exact_phase_values() -> np.ndarray:
    """Return all 6^9 determinants, using integer operations only."""

    shape = (6,) * len(CYCLE_EDGES)
    real = np.zeros(shape, dtype=np.int64)
    omega = np.zeros(shape, dtype=np.int64)
    for exponents, coefficient in determinant_polynomial().items():
        real[exponents] = coefficient
    for axis in range(len(CYCLE_EDGES)):
        real, omega = _transform_axis(real, omega, axis)
        assert int(np.max(np.abs(real))) < 2**62
        assert int(np.max(np.abs(omega))) < 2**62
    # Hermitian determinants are rational integers.
    assert not np.any(omega)
    return real


def certificate() -> dict[str, object]:
    values = exact_phase_values()
    internal_maximum = int(np.max(values)) * 3**ORDER
    assert internal_maximum == 7464960
    full_maximum = internal_maximum * 15**9
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    assert full_maximum < q153
    return {
        "support": "K6 minus one edge plus nine isolates",
        "switching_normalized_phase_assignments": 6**len(CYCLE_EDGES),
        "distinct_scaled_internal_determinants": int(np.unique(values).size),
        "scaled_internal_maximum": int(np.max(values)),
        "internal_maximum": internal_maximum,
        "full_maximum": full_maximum,
        "over_q153": float(full_maximum / q153),
        "theorem": "Every Q=126 Gram with this support lies below the Q=153 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
