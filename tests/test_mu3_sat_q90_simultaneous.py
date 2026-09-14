from order15_mu3.scripts.sat_q90_simultaneous import build, target_norm


def test_q90_simultaneous_model_builds() -> None:
    cnf, pool, _x = build()
    assert target_norm(0, 1) == 9
    assert target_norm(0, 2) == 0
    assert target_norm(0, 13) == 3
    assert pool.top > 0 and cnf.clauses
