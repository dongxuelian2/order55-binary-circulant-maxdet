from maxdet.mu3 import Eisenstein, determinant
from order15_mu3.scripts.sat_q90_top import TARGET, build, target_gram


def test_q90_top_target_and_model() -> None:
    gram = target_gram()
    assert determinant(gram) == Eisenstein(TARGET)
    cnf, pool, _x = build()
    assert pool.top > 0 and cnf.clauses
