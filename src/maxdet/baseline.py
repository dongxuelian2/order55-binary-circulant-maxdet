"""Minimal random-search scaffold for smoke testing the experiment pipeline."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .core import exact_determinant, fast_determinant, random_pm1_matrix


@dataclass(frozen=True)
class BaselineResult:
    """Best candidate returned by a deterministic random baseline run."""

    n: int
    samples: int
    seed: int
    best_fast_score: float
    exact_determinant: int
    matrix: np.ndarray

    def as_dict(self) -> dict[str, object]:
        """Return JSON-compatible result fields without run metadata."""

        return {
            "n": self.n,
            "samples": self.samples,
            "seed": self.seed,
            "best_fast_score": self.best_fast_score,
            "exact_determinant": self.exact_determinant,
            "matrix": self.matrix.tolist(),
        }


def run_random_baseline(n: int, samples: int, seed: int) -> BaselineResult:
    """Generate, score, and exactly certify a small random sample.

    The NumPy determinant is used only for ranking. The returned candidate is
    certified with SymPy before the result is returned.
    """

    if not isinstance(n, (int, np.integer)) or isinstance(n, (bool, np.bool_)):
        raise TypeError("n must be an integer")
    if n < 1:
        raise ValueError("n must be positive")
    if not isinstance(samples, (int, np.integer)) or isinstance(
        samples, (bool, np.bool_)
    ):
        raise TypeError("samples must be an integer")
    if samples < 1:
        raise ValueError("samples must be positive")
    if not isinstance(seed, (int, np.integer)) or isinstance(seed, (bool, np.bool_)):
        raise TypeError("seed must be an integer")

    generator = np.random.default_rng(int(seed))
    best_matrix: np.ndarray | None = None
    best_score = float("-inf")

    for _ in range(int(samples)):
        candidate = random_pm1_matrix(int(n), generator)
        score = abs(fast_determinant(candidate))
        if score > best_score:
            best_score = score
            best_matrix = candidate

    assert best_matrix is not None
    certified = exact_determinant(best_matrix)
    return BaselineResult(
        n=int(n),
        samples=int(samples),
        seed=int(seed),
        best_fast_score=float(best_score),
        exact_determinant=certified,
        matrix=best_matrix.copy(),
    )
