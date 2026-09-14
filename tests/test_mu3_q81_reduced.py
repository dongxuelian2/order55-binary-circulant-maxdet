from order15_mu3.scripts.sat_q81_reduced import build
from pysat.solvers import Solver


def test_q81_reduced_block_model_builds() -> None:
    cnf, pool, cell = build()
    assert cnf.clauses and pool.top > 0
    assert cell(0, 0, 0) == cell(2, 2, 0)
    assert cell(3, 3, 0) == cell(4, 4, 0)
    assert cell(3, 4, 0) == cell(4, 3, 0)


def test_q81_reduced_block_model_is_unsatisfiable() -> None:
    cnf, _pool, _cell = build()
    with Solver(name="glucose42", bootstrap_with=cnf) as solver:
        assert solver.solve() is False
