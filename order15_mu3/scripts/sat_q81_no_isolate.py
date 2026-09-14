"""Exact SAT test for the final no-isolate all-same-color Q=81 profile.

After independent row and column phase gauges, both Gram supports must be six
disjoint +3 edges and one -3 triangle.  This model enforces those complete row
and column Gram matrices simultaneously.
"""

from __future__ import annotations

import argparse
import json
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import gram_from_exponents


ORDER = 15
PAIR_EDGES = {(0, 1), (2, 3), (4, 5), (6, 7), (8, 9), (10, 11)}
TRIANGLE_EDGES = {(12, 13), (12, 14), (13, 14)}


def target_counts(first: int, second: int) -> tuple[int, int, int]:
    edge = (first, second)
    if edge in PAIR_EDGES:
        return 7, 4, 4  # Gram entry +3.
    if edge in TRIANGLE_EDGES:
        return 3, 6, 6  # Gram entry -3.
    return 5, 5, 5


def build():
    pool = IDPool()
    cnf = CNF()

    def x(row: int, column: int, phase: int) -> int:
        return pool.id(("x", row, column, phase))

    for row in range(ORDER):
        for column in range(ORDER):
            cnf.extend(CardEnc.equals([x(row, column, phase) for phase in range(3)], 1, vpool=pool, encoding=EncType.pairwise))

    # Canonicalizing both Gram matrices has already consumed the relative row
    # and column phases.  Only the overall scalar phase remains universally
    # available, so fix one entry and no more.
    cnf.append([x(0, 0, 0)])

    def constrain_pair(axis: str, first: int, second: int) -> None:
        differences = [[], [], []]
        for position in range(ORDER):
            for difference in range(3):
                y = pool.id((axis, first, second, position, difference))
                differences[difference].append(y)
                for phase in range(3):
                    if axis == "row_difference":
                        left = x(first, position, phase)
                        right = x(second, position, (phase - difference) % 3)
                    else:
                        left = x(position, first, phase)
                        right = x(position, second, (phase - difference) % 3)
                    cnf.append([-left, -right, y])
                    cnf.append([-y, -left, right])
        for literals, count in zip(differences, target_counts(first, second)):
            cnf.extend(CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter))

    for first in range(ORDER):
        for second in range(first + 1, ORDER):
            constrain_pair("row_difference", first, second)
            constrain_pair("column_difference", first, second)
    return cnf, pool, x


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", default="glucose42")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    args = parser.parse_args()
    cnf, pool, x = build()
    with Solver(name=args.solver, bootstrap_with=cnf) as solver:
        timer = Timer(args.timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "model": "Q=81 no-isolate K3+6K2 row and column Gram matrices",
            "solver": args.solver,
            "timeout_seconds": args.timeout_seconds,
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
            columns = [[rows[row][column] for row in range(ORDER)] for column in range(ORDER)]
            column_gram = gram_from_exponents(columns)
            for gram in (row_gram, column_gram):
                for first in range(ORDER):
                    for second in range(first + 1, ORDER):
                        expected = 9 if (first, second) in PAIR_EDGES | TRIANGLE_EDGES else 0
                        assert gram[first][second].norm() == expected
            payload["rows"] = rows
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
