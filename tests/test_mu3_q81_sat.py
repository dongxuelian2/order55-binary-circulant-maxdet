from order15_mu3.scripts.sat_q81_no_isolate import PAIR_EDGES, TRIANGLE_EDGES, build, target_counts


def test_q81_no_isolate_complete_gram_model_builds() -> None:
    assert len(PAIR_EDGES) == 6
    assert len(TRIANGLE_EDGES) == 3
    assert target_counts(0, 1) == (7, 4, 4)
    assert target_counts(12, 13) == (3, 6, 6)
    assert target_counts(0, 2) == (5, 5, 5)
    cnf, pool, _x = build()
    assert cnf.clauses and pool.top > 0
