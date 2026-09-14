"""Exact simultaneous row/column support SAT model for the Q=90 size-13 case."""

from __future__ import annotations

import argparse
import json
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import determinant, gram_from_exponents


ORDER = 15
LARGE = range(13)
HUBS = (13, 14)
NORM_COUNTS = {
    0: ((5, 5, 5),),
    3: ((4, 5, 6), (4, 6, 5), (5, 4, 6), (5, 6, 4), (6, 4, 5), (6, 5, 4)),
    9: ((7, 4, 4), (4, 7, 4), (4, 4, 7), (6, 6, 3), (6, 3, 6), (3, 6, 6)),
}


def target_norm(first: int, second: int) -> int:
    if first < 13 and second < 13:
        return 9 if (first, second) == (0, 1) else 0
    return 3


def _guarded_equals(cnf: CNF, pool: IDPool, literals: list[int], count: int, guard: int) -> None:
    encoded = CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter)
    cnf.extend([[-guard, *clause] for clause in encoded.clauses])


def _axis_lex(cnf: CNF, pool: IDPool, x, axis: str, first: int, second: int) -> None:
    prefix = pool.id((axis, "lex-prefix", first, second, 0))
    cnf.append([prefix])
    for position in range(ORDER):
        same = pool.id((axis, "lex-same", first, second, position))
        left = lambda phase: x(first, position, phase) if axis == "row" else x(position, first, phase)
        right = lambda phase: x(second, position, phase) if axis == "row" else x(position, second, phase)
        for phase in range(3):
            cnf.append([-left(phase), -right(phase), same])
            cnf.append([-same, -left(phase), right(phase)])
        for left_phase in range(3):
            for right_phase in range(left_phase):
                cnf.append([-prefix, -left(left_phase), -right(right_phase)])
        following = pool.id((axis, "lex-prefix", first, second, position + 1))
        cnf.append([-following, prefix])
        cnf.append([-following, same])
        cnf.append([-prefix, -same, following])
        prefix = following


def build():
    pool = IDPool()
    cnf = CNF()

    def x(row: int, column: int, phase: int) -> int:
        return pool.id(("x", row, column, phase))

    for row in range(ORDER):
        for column in range(ORDER):
            cnf.extend(CardEnc.equals(
                [x(row, column, phase) for phase in range(3)],
                1, vpool=pool, encoding=EncType.pairwise,
            ))
    # Independent row and column phase switching dephases the first row and
    # first column without changing either Gram support.
    for position in range(ORDER):
        cnf.append([x(0, position, 0)])
        cnf.append([x(position, 0, 0)])

    for axis in ("row", "column"):
        for first in range(ORDER):
            for second in range(first + 1, ORDER):
                differences = [[], [], []]
                for position in range(ORDER):
                    for difference in range(3):
                        indicator = pool.id((axis, "difference", first, second, position, difference))
                        differences[difference].append(indicator)
                        for phase in range(3):
                            if axis == "row":
                                left = x(first, position, phase)
                                right = x(second, position, (phase - difference) % 3)
                            else:
                                left = x(position, first, phase)
                                right = x(position, second, (phase - difference) % 3)
                            cnf.append([-left, -right, indicator])
                            cnf.append([-indicator, -left, right])
                patterns = NORM_COUNTS[target_norm(first, second)]
                if len(patterns) == 1:
                    for literals, count in zip(differences, patterns[0]):
                        cnf.extend(CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter))
                else:
                    selectors = [pool.id((axis, "pattern", first, second, index)) for index in range(len(patterns))]
                    cnf.extend(CardEnc.equals(selectors, 1, vpool=pool, encoding=EncType.pairwise))
                    for selector, counts in zip(selectors, patterns):
                        for literals, count in zip(differences, counts):
                            _guarded_equals(cnf, pool, literals, count, selector)

    # The eleven isolated large-class support vertices are interchangeable
    # independently on the row and column axes.
    for axis in ("row", "column"):
        for first in range(2, 12):
            _axis_lex(cnf, pool, x, axis, first, first + 1)
    return cnf, pool, x


def solve(solver_name: str = "cadical195", timeout_seconds: float = 300) -> dict[str, object]:
    cnf, pool, x = build()
    with Solver(name=solver_name, bootstrap_with=cnf) as solver:
        timer = Timer(timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "model": "simultaneous Q=90 (13,1,1) row/column Gram norm supports",
            "solver": solver_name,
            "variables": pool.top,
            "clauses": len(cnf.clauses),
            "result": "sat" if result is True else "unsat" if result is False else "unknown_timeout",
        }
        if result is True:
            model = {literal for literal in solver.get_model() if literal > 0}
            rows = [
                [next(phase for phase in range(3) if x(row, column, phase) in model) for column in range(ORDER)]
                for row in range(ORDER)
            ]
            row_gram = gram_from_exponents(rows)
            column_gram = gram_from_exponents([list(column) for column in zip(*rows)])
            assert all(
                gram[first][second].norm() == target_norm(first, second)
                for gram in (row_gram, column_gram)
                for first in range(ORDER) for second in range(first + 1, ORDER)
            )
            value = determinant(row_gram)
            assert value == determinant(column_gram)
            payload["determinant"] = value.a
            payload["rows"] = rows
        return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    args = parser.parse_args()
    print(json.dumps(solve(args.solver, args.timeout_seconds), indent=2))


if __name__ == "__main__":
    main()
