from __future__ import annotations

import pytest

from maxdet.core import fast_logabs_determinant
from maxdet.twocirculant import (
    bordered_two_circulant,
    parse_sign_word,
    spectral_logabsdet,
)


@pytest.mark.parametrize("matrix_type", [1, 3])
def test_spectral_score_matches_direct_slogdet(matrix_type: int) -> None:
    a = parse_sign_word("+-+-+")
    b = parse_sign_word("++---")
    matrix = bordered_two_circulant(a, b, matrix_type)
    assert matrix.shape == (11, 11)
    assert spectral_logabsdet(a, b, matrix_type) == pytest.approx(
        fast_logabs_determinant(matrix)
    )
