from order15_mu3.scripts.sat_q78_paw_simultaneous import build


def test_simultaneous_q78_paw_encoding_shape() -> None:
    cnf, pool, _ = build(0, 0)
    assert pool.top > 10000
    assert len(cnf.clauses) > 100000
