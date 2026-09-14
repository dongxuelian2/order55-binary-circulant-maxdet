"""Vectorized scan of all 3^14 dephased order-15 circulant first rows.

The FFT stage is a discovery/ranking pass. Finalists are recomputed exactly in
Z[omega], so the reported construction values are certified; floating-point
ranking alone is not claimed as a proof of circulant optimality.
"""

from __future__ import annotations

import argparse
import json

import numpy as np

from maxdet.mu3 import determinant, exponent_matrix
from maxdet.mu3_search import COMPLEX_ROOTS


def digits_base_three(start: int, count: int) -> np.ndarray:
    values = np.arange(start, start + count, dtype=np.int64)
    digits = np.zeros((count, 15), dtype=np.int8)
    for column in range(1, 15):
        digits[:, column] = values % 3
        values //= 3
    return digits


def circulant(first_row: np.ndarray) -> list[list[int]]:
    return [np.roll(first_row, shift).tolist() for shift in range(len(first_row))]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--batch-size", type=int, default=100_000)
    parser.add_argument("--finalists", type=int, default=128)
    args = parser.parse_args()
    total = 3**14
    retained: list[tuple[float, int, list[int]]] = []
    for start in range(0, total, args.batch_size):
        count = min(args.batch_size, total - start)
        exponents = digits_base_three(start, count)
        rows = COMPLEX_ROOTS[exponents]
        eigenvalues = np.fft.fft(rows, axis=1)
        with np.errstate(divide="ignore"):
            scores = np.sum(np.log(np.abs(eigenvalues)), axis=1)
        keep = min(args.finalists, count)
        indices = np.argpartition(scores, -keep)[-keep:]
        retained.extend((float(scores[index]), start + int(index), exponents[index].tolist()) for index in indices)
        retained = sorted(retained, reverse=True)[: args.finalists]

    exact: list[tuple[int, int, int, int, list[int]]] = []
    for _score, code, first_row in retained:
        value = determinant(exponent_matrix(circulant(np.asarray(first_row, dtype=np.int8))))
        exact.append((value.norm(), value.a, value.b, code, first_row))
    exact.sort(reverse=True)
    norm, a, b, code, first_row = exact[0]
    print(json.dumps({
        "scanned": total,
        "floating_ranked_finalists": len(exact),
        "best_exact_norm_among_finalists": norm,
        "determinant": [a, b],
        "base3_code": code,
        "first_row": first_row,
        "proof_status": "exact value for finalist; floating scan is not an exact optimality certificate",
    }, indent=2))


if __name__ == "__main__":
    main()
