"""Exactly extend a certified nine-row orthogonal seed by one more row."""

from __future__ import annotations

import json

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import gram_from_exponents


SEED = (
    (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
    (0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2),
    (1, 1, 2, 2, 2, 0, 0, 1, 1, 1, 0, 0, 0, 2, 2),
    (0, 2, 2, 1, 1, 1, 0, 1, 0, 2, 1, 2, 0, 2, 0),
    (0, 2, 1, 2, 1, 1, 1, 0, 2, 0, 0, 0, 1, 2, 2),
    (0, 2, 0, 1, 2, 2, 1, 2, 0, 1, 0, 1, 0, 1, 2),
    (0, 1, 1, 0, 2, 0, 1, 2, 0, 2, 1, 0, 2, 2, 1),
    (0, 2, 1, 2, 1, 0, 0, 2, 1, 1, 2, 1, 2, 0, 0),
    (0, 1, 2, 2, 0, 1, 2, 2, 0, 1, 0, 2, 1, 0, 1),
)


def main() -> None:
    gram = gram_from_exponents(SEED)
    assert all(gram[i][j].norm() == 0 for i in range(len(SEED)) for j in range(i))
    pool = IDPool()
    cnf = CNF()

    def x(column: int, phase: int) -> int:
        return pool.id((column, phase))

    for column in range(15):
        cnf.extend(CardEnc.equals([x(column, phase) for phase in range(3)], 1, vpool=pool, encoding=EncType.pairwise))
    cnf.append([x(0, 0)])
    for seed_row in SEED:
        for difference in range(3):
            literals = [x(column, (seed_row[column] + difference) % 3) for column in range(15)]
            cnf.extend(CardEnc.equals(literals, 5, vpool=pool, encoding=EncType.seqcounter))

    with Solver(name="cadical195", bootstrap_with=cnf) as solver:
        result = solver.solve()
        payload = {
            "seed_rows": len(SEED),
            "result": "sat" if result else "unsat",
            "variables": pool.top,
            "clauses": len(cnf.clauses),
        }
        if result:
            model = {literal for literal in solver.get_model() if literal > 0}
            extension = [next(phase for phase in range(3) if x(column, phase) in model) for column in range(15)]
            rows = [list(row) for row in SEED] + [extension]
            gram = gram_from_exponents(rows)
            assert all(gram[i][j].norm() == 0 for i in range(len(rows)) for j in range(i))
            payload["extension"] = extension
            payload["rows"] = rows
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
