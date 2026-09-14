import numpy as np

from order15_mu3.scripts.staged_friendship_search import norms_to


def test_vectorized_gram_norms() -> None:
    rows = np.asarray([[0, 1, 2], [0, 0, 1], [2, 2, 0]], dtype=np.int8)
    target = np.asarray([0, 1, 2], dtype=np.int8)
    assert norms_to(rows, target).tolist() == [9, 3, 3]
