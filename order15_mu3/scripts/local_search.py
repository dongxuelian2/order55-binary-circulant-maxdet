"""Iterated single-entry search with exact certification of promising optima."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from maxdet.mu3 import dephase, determinant, exponent_matrix, validate_exponents
from maxdet.mu3_search import (
    alternating_row_column_ascent,
    perturb,
    random_dephased,
    single_entry_ascent,
)


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "data" / "literature_matrices" / "nunez_ponasso_m15.json"


def exact_norm(x: np.ndarray) -> tuple[int, int, int]:
    value = determinant(exponent_matrix(x.tolist()))
    return value.norm(), value.a, value.b


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--iterations", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1503)
    parser.add_argument("--random-every", type=int, default=20)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    published = validate_exponents(json.loads(SOURCE.read_text(encoding="utf-8"))["matrix"], 15)
    incumbent = np.asarray(dephase(published), dtype=np.int8)
    best_norm, best_a, best_b = exact_norm(incumbent)
    best_log = 0.5 * np.log(best_norm)
    basins = 0
    singular = 0
    for iteration in range(args.iterations):
        if iteration and iteration % args.random_every == 0:
            start = random_dephased(15, rng)
        else:
            start = perturb(incumbent, int(rng.integers(2, 13)), rng)
        try:
            result = single_entry_ascent(start)
            result = alternating_row_column_ascent(result.exponents)
        except np.linalg.LinAlgError:
            singular += 1
            continue
        basins += 1
        if result.logabsdet > best_log + 1e-8:
            norm, a, b = exact_norm(result.exponents)
            if norm > best_norm:
                incumbent = result.exponents.copy()
                best_norm, best_a, best_b = norm, a, b
                best_log = 0.5 * np.log(best_norm)
                print(f"record iteration={iteration} norm={best_norm} det=({best_a},{best_b})", flush=True)
    output = {
        "seed": args.seed,
        "iterations": args.iterations,
        "completed_basins": basins,
        "singular_starts": singular,
        "best_norm": best_norm,
        "determinant": [best_a, best_b],
        "matrix": incumbent.tolist(),
    }
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
