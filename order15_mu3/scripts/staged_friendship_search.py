"""Staged exact-row augmentation search for the all-norm-3 F7 support."""

from __future__ import annotations

import argparse
import json
import time
from itertools import combinations

import numpy as np

from maxdet.mu3 import determinant, exponent_matrix, gram_from_exponents


def canonical_row_universe() -> np.ndarray:
    """Generate all 630630 rows with exponent counts (4,6,5)."""

    total = 630_630
    rows = np.full((total, 15), 2, dtype=np.int8)
    index = 0
    positions = tuple(range(15))
    for zeros in combinations(positions, 4):
        zero_set = set(zeros)
        remaining = tuple(position for position in positions if position not in zero_set)
        for ones in combinations(remaining, 6):
            rows[index, list(zeros)] = 0
            rows[index, list(ones)] = 1
            index += 1
    assert index == total
    return rows


def norms_to(rows: np.ndarray, target: np.ndarray) -> np.ndarray:
    differences = (rows - target) % 3
    c0 = np.sum(differences == 0, axis=1)
    c1 = np.sum(differences == 1, axis=1)
    c2 = rows.shape[1] - c0 - c1
    a = c0 - c2
    b = c1 - c2
    return a * a - a * b + b * b


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--seconds", type=float, default=180.0)
    parser.add_argument("--seed", type=int, default=9157)
    parser.add_argument("--first-samples", type=int, default=24)
    parser.add_argument("--mate-samples", type=int, default=12)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    universe = canonical_row_universe()
    first = np.asarray((0,) * 4 + (1,) * 6 + (2,) * 5, dtype=np.int8)
    relation_first = norms_to(universe, first)
    possible_mates = universe[relation_first == 3]
    if len(possible_mates) == 0:
        print(json.dumps({
            "universe_size": len(universe),
            "first_leaf_possible_mates": 0,
            "result": "unsat",
            "proof_status": "exact finite enumeration; explained by the row-sum color congruence",
        }, indent=2))
        return
    deadline = time.monotonic() + args.seconds
    trials = 0
    deepest = 1
    best_pool_sizes: list[int] = []
    solution: list[np.ndarray] | None = None
    while time.monotonic() < deadline and solution is None:
        trials += 1
        second = possible_mates[int(rng.integers(len(possible_mates)))]
        pool = universe[(relation_first == 0) & (norms_to(universe, second) == 0)]
        chosen = [first, second]
        pool_sizes = [len(pool)]
        while len(chosen) < 14 and len(pool):
            sample_count = min(args.first_samples, len(pool))
            first_indices = rng.choice(len(pool), size=sample_count, replace=False)
            best: tuple[int, np.ndarray, np.ndarray, np.ndarray] | None = None
            for first_index in first_indices:
                candidate = pool[int(first_index)]
                relation_candidate = norms_to(pool, candidate)
                mate_indices = np.flatnonzero(relation_candidate == 3)
                if not len(mate_indices):
                    continue
                if len(mate_indices) > args.mate_samples:
                    mate_indices = rng.choice(mate_indices, size=args.mate_samples, replace=False)
                orth_candidate = relation_candidate == 0
                for mate_index in mate_indices:
                    mate = pool[int(mate_index)]
                    next_pool = pool[orth_candidate & (norms_to(pool, mate) == 0)]
                    score = len(next_pool)
                    if best is None or score > best[0]:
                        best = (score, candidate.copy(), mate.copy(), next_pool)
            if best is None:
                break
            _, candidate, mate, pool = best
            chosen.extend((candidate, mate))
            pool_sizes.append(len(pool))
            if len(chosen) // 2 > deepest:
                deepest = len(chosen) // 2
                best_pool_sizes = pool_sizes[:]
                print(f"deepest_pairs={deepest} trial={trials} pools={best_pool_sizes}", flush=True)
        if len(chosen) == 14:
            solution = chosen

    payload: dict[str, object] = {
        "seed": args.seed,
        "seconds": args.seconds,
        "universe_size": len(universe),
        "first_leaf_possible_mates": len(possible_mates),
        "trials": trials,
        "deepest_pairs": deepest,
        "best_pool_sizes": best_pool_sizes,
        "result": "found" if solution is not None else "not_found",
        "proof_status": "a found matrix is exact; not_found is heuristic",
    }
    if solution is not None:
        matrix = np.vstack((np.zeros(15, dtype=np.int8), np.asarray(solution)))
        gram = gram_from_exponents(matrix.tolist())
        for i in range(1, 15):
            assert gram[0][i].norm() == 3
        for i in range(1, 15):
            for j in range(i + 1, 15):
                assert gram[i][j].norm() == (3 if (i - 1) // 2 == (j - 1) // 2 else 0)
        value = determinant(exponent_matrix(matrix.tolist()))
        payload.update({"determinant": [value.a, value.b], "determinant_norm": value.norm(), "matrix": matrix.tolist()})
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
