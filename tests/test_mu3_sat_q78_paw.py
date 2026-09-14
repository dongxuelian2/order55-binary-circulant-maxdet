from maxdet.mu3 import Eisenstein

from order15_mu3.scripts.sat_q78_paw import LABELINGS, MAXIMUM, PAW_EDGES, build, target_gram


def test_q78_paw_encoding_shape() -> None:
    assert len(LABELINGS) == 4
    assert len(set(LABELINGS)) == 4
    assert MAXIMUM == 308697676962890625
    assert all(target_gram(*edge, 0).norm() == 9 for edge in PAW_EDGES)
    assert target_gram(4, 14, 0) == Eisenstein(1, -1)
    cnf, pool, _ = build(0)
    assert pool.top > 0 and len(cnf.clauses) > 0
