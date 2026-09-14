"""Simultaneous row/column Gram SAT test at the Q=78 determinant cap."""

from __future__ import annotations

import argparse
import json
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import determinant, gram_from_exponents

try:
    from order15_mu3.scripts.sat_q78_paw import LABELINGS, MAXIMUM, ORDER, count_lookup, target_gram
except ModuleNotFoundError:
    from sat_q78_paw import LABELINGS, MAXIMUM, ORDER, count_lookup, target_gram


def add_axis_lex(cnf: CNF, pool: IDPool, x, axis: str, first: int, second: int) -> None:
    prefix = pool.id((axis, "lex_prefix", first, second, 0))
    cnf.append([prefix])
    for position in range(ORDER):
        same = pool.id((axis, "lex_same", first, second, position))
        left_vars = lambda phase: x(first, position, phase) if axis == "row" else x(position, first, phase)
        right_vars = lambda phase: x(second, position, phase) if axis == "row" else x(position, second, phase)
        for phase in range(3):
            cnf.append([-left_vars(phase), -right_vars(phase), same])
            cnf.append([-same, -left_vars(phase), right_vars(phase)])
        for left_phase in range(3):
            for right_phase in range(left_phase):
                cnf.append([-prefix, -left_vars(left_phase), -right_vars(right_phase)])
        next_prefix = pool.id((axis, "lex_prefix", first, second, position + 1))
        cnf.append([-next_prefix, prefix])
        cnf.append([-next_prefix, same])
        cnf.append([-prefix, -same, next_prefix])
        prefix = next_prefix


def build(row_labeling: int, column_labeling: int):
    counts = count_lookup()
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
    cnf.append([x(0, 0, 0)])

    for axis, labeling in (("row", row_labeling), ("column", column_labeling)):
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
                target = counts[target_gram(first, second, labeling)]
                for literals, count in zip(differences, target):
                    cnf.extend(CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter))

    # The ten zero-degree paw vertices are interchangeable on both axes.
    for axis in ("row", "column"):
        for first in range(4, 13):
            add_axis_lex(cnf, pool, x, axis, first, first + 1)
    return cnf, pool, x


def solve(row_labeling: int, column_labeling: int, solver_name: str, timeout_seconds: float):
    cnf, pool, x = build(row_labeling, column_labeling)
    with Solver(name=solver_name, bootstrap_with=cnf) as solver:
        timer = Timer(timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "model": "simultaneous Q=78 maximum paw row and column Grams",
            "row_labeling": row_labeling,
            "column_labeling": column_labeling,
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
            assert determinant(row_gram).a == MAXIMUM
            assert determinant(column_gram).a == MAXIMUM
            payload["rows"] = rows
        return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("row_labeling", type=int, choices=range(len(LABELINGS)))
    parser.add_argument("column_labeling", type=int, choices=range(len(LABELINGS)))
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    args = parser.parse_args()
    print(json.dumps(solve(args.row_labeling, args.column_labeling, args.solver, args.timeout_seconds), indent=2))


if __name__ == "__main__":
    main()
