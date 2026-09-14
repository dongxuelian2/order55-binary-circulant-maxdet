"""Exact verification of the published order-15 benchmark matrix."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path

from maxdet.mu3 import determinant, exponent_matrix, gram_from_exponents, validate_exponents


ROOT = Path(__file__).resolve().parents[1]
MATRIX_PATH = ROOT / "data" / "literature_matrices" / "nunez_ponasso_m15.json"


def main() -> None:
    record = json.loads(MATRIX_PATH.read_text(encoding="utf-8"))
    exponents = validate_exponents(record["matrix"], order=15)
    det_h = determinant(exponent_matrix(exponents))
    gram = gram_from_exponents(exponents)
    det_gram = determinant(gram)
    expected_norm = 2**22 * 3**20 * 19
    assert det_h.norm() == expected_norm
    assert det_gram.a == expected_norm and det_gram.b == 0

    off_diagonal_norms = Counter(
        gram[i][j].norm()
        for i in range(15)
        for j in range(i)
        if gram[i][j].norm() != 0
    )
    zero_pairs = sum(gram[i][j].norm() == 0 for i in range(15) for j in range(i))
    payload = {
        "determinant": {"a": det_h.a, "b": det_h.b},
        "squared_absolute_determinant": expected_norm,
        "factorization": "2^22 * 3^20 * 19",
        "normalized_squared_determinant": expected_norm // 3**14,
        "normalized_factorization": "2^22 * 3^6 * 19",
        "zero_off_diagonal_pairs": zero_pairs,
        "nonzero_off_diagonal_norm_histogram": dict(sorted(off_diagonal_norms.items())),
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
