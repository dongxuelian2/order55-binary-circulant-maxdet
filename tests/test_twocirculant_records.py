from __future__ import annotations

import numpy as np
import pytest

from maxdet.core import fast_logabs_determinant
from maxdet.twocirculant import bordered_two_circulant, parse_sign_word, spectral_logabsdet


@pytest.mark.parametrize(
    "a,b,matrix_type",
    [
        (
            "-+-++-++++-+-+--++++--++-+-+++-+--------+++-+---++-",
            "-+++++--+--+-+----+++-+-++++-+--++-+++---++-++--+++",
            1,
        ),
        (
            "++--++--+---++++---++---+++-++--+-+--++-++++-+-+++++-+-",
            "----++++++++-+-++-++-+-+-++---+-++-----+-+++-+--+++----",
            3,
        ),
    ],
)
def test_public_record_spectral_score(a: str, b: str, matrix_type: int) -> None:
    first, second = parse_sign_word(a), parse_sign_word(b)
    matrix = bordered_two_circulant(first, second, matrix_type)
    assert spectral_logabsdet(first, second, matrix_type) == pytest.approx(
        fast_logabs_determinant(matrix), abs=1e-10
    )
    assert np.all((matrix == -1) | (matrix == 1))
