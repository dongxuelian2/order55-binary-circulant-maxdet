"""Constraint-preserving search for an all-norm-3 F7 support realization."""

from __future__ import annotations

import argparse
import json
import math
import time

import numpy as np

from maxdet.mu3 import determinant, exponent_matrix, gram_from_exponents


ORDER = 15
PAIRS = {frozenset((1 + 2 * k, 2 + 2 * k)) for k in range(7)}
CANONICAL_COUNTS = (4, 6, 5)


def gram_norm(first: np.ndarray, second: np.ndarray) -> int:
    counts = np.bincount((first - second) % 3, minlength=3)
    a = int(counts[0] - counts[2])
    b = int(counts[1] - counts[2])
    return a * a - a * b + b * b


def pair_penalty(first: int, second: int, rows: np.ndarray) -> int:
    target = 3 if frozenset((first, second)) in PAIRS else 0
    return abs(gram_norm(rows[first], rows[second]) - target)


def total_penalty(rows: np.ndarray) -> int:
    return sum(pair_penalty(i, j, rows) for i in range(1, ORDER) for j in range(i + 1, ORDER))


def random_state(rng: np.random.Generator) -> np.ndarray:
    rows = np.zeros((ORDER, ORDER), dtype=np.int8)
    rows[1] = np.asarray((0,) * 4 + (1,) * 6 + (2,) * 5, dtype=np.int8)
    for row in range(2, ORDER):
        values = np.repeat(np.arange(3, dtype=np.int8), CANONICAL_COUNTS)
        rng.shuffle(values)
        rows[row] = values
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=120.0)
    parser.add_argument("--seed", type=int, default=7315)
    parser.add_argument("--moves-per-restart", type=int, default=200_000)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    deadline = time.monotonic() + args.seconds
    best_penalty = 10**9
    best_rows: np.ndarray | None = None
    moves = restarts = 0
    while time.monotonic() < deadline and best_penalty:
        rows = random_state(rng)
        penalty = total_penalty(rows)
        restarts += 1
        for step in range(args.moves_per_restart):
            if time.monotonic() >= deadline or penalty == 0:
                break
            row = int(rng.integers(2, ORDER))
            first, second = rng.choice(ORDER, size=2, replace=False)
            if rows[row, first] == rows[row, second]:
                continue
            old_local = sum(pair_penalty(min(row, other), max(row, other), rows) for other in range(1, ORDER) if other != row)
            rows[row, first], rows[row, second] = rows[row, second], rows[row, first]
            new_local = sum(pair_penalty(min(row, other), max(row, other), rows) for other in range(1, ORDER) if other != row)
            delta = new_local - old_local
            temperature = max(0.05, 12.0 * (1.0 - step / args.moves_per_restart))
            if delta <= 0 or rng.random() < math.exp(-delta / temperature):
                penalty += delta
            else:
                rows[row, first], rows[row, second] = rows[row, second], rows[row, first]
            moves += 1
            if penalty < best_penalty:
                best_penalty = penalty
                best_rows = rows.copy()
                print(f"best penalty={best_penalty} restarts={restarts} moves={moves}", flush=True)
    assert best_rows is not None
    payload: dict[str, object] = {
        "seed": args.seed,
        "seconds": args.seconds,
        "restarts": restarts,
        "moves": moves,
        "best_penalty": best_penalty,
        "proof_status": "heuristic unless best_penalty is zero; zero is exactly rechecked",
        "leaf_pair_norm_histogram": {
            str(norm): sum(
                gram_norm(best_rows[i], best_rows[j]) == norm
                for i in range(1, ORDER)
                for j in range(i + 1, ORDER)
            )
            for norm in sorted({gram_norm(best_rows[i], best_rows[j]) for i in range(1, ORDER) for j in range(i + 1, ORDER)})
        },
        "violating_pairs": [
            [i, j, gram_norm(best_rows[i], best_rows[j]), 3 if frozenset((i, j)) in PAIRS else 0]
            for i in range(1, ORDER)
            for j in range(i + 1, ORDER)
            if gram_norm(best_rows[i], best_rows[j]) != (3 if frozenset((i, j)) in PAIRS else 0)
        ],
        "best_matrix": best_rows.tolist(),
    }
    if best_penalty == 0:
        gram = gram_from_exponents(best_rows.tolist())
        exact_penalty = sum(abs(gram[i][j].norm() - (3 if frozenset((i, j)) in PAIRS else 0)) for i in range(1, ORDER) for j in range(i + 1, ORDER))
        assert exact_penalty == 0
        value = determinant(exponent_matrix(best_rows.tolist()))
        payload.update({"determinant": [value.a, value.b], "determinant_norm": value.norm(), "matrix": best_rows.tolist()})
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
