"""Small, exactness-aware primitives for ±1 determinant experiments."""

from .baseline import BaselineResult, run_random_baseline
from .core import (
    exact_determinant,
    fast_determinant,
    fast_logabs_determinant,
    negate_column,
    negate_row,
    normalize_pm1_matrix,
    random_pm1_matrix,
    swap_columns,
    swap_rows,
    sylvester_hadamard,
    validate_pm1_matrix,
)

__all__ = [
    "BaselineResult",
    "exact_determinant",
    "fast_determinant",
    "fast_logabs_determinant",
    "negate_column",
    "negate_row",
    "normalize_pm1_matrix",
    "random_pm1_matrix",
    "run_random_baseline",
    "swap_columns",
    "swap_rows",
    "sylvester_hadamard",
    "validate_pm1_matrix",
]

__version__ = "0.1.0"
