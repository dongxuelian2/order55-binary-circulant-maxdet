"""Reduced exact SAT model for the last no-isolate Q=81 profile.

The -1 triangle phase and A H = H B force a 3+6*2 block form.  The triangle
corner is all ones, each outer row/column pair repeats a permutation of the
three roots on the triangle, and every outer 2x2 block is [[x,y],[y,x]].
This model enforces the complete row and column Gram counts on that form.
"""

from __future__ import annotations

import argparse
import json
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

try:
    from order15_mu3.scripts.sat_q81_no_isolate import ORDER, target_counts
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from sat_q81_no_isolate import ORDER, target_counts


def build():
    pool = IDPool()
    cnf = CNF()

    def symbol(kind: str, first: int, second: int, phase: int) -> int:
        return pool.id((kind, first, second, phase))

    constant = [pool.id(("constant", phase)) for phase in range(3)]
    cnf.append([constant[0]])
    cnf.append([-constant[1]])
    cnf.append([-constant[2]])

    # r[a,j] and q[b,i] are permutations of 0,1,2.  x[a,b],y[a,b]
    # are unrestricted third-root exponents.
    for kind in ("r", "q"):
        for block in range(6):
            for position in range(3):
                cnf.extend(CardEnc.equals([symbol(kind, block, position, phase) for phase in range(3)], 1, vpool=pool, encoding=EncType.pairwise))
            for phase in range(3):
                cnf.extend(CardEnc.equals([symbol(kind, block, position, phase) for position in range(3)], 1, vpool=pool, encoding=EncType.pairwise))
    for kind in ("x", "y"):
        for row_block in range(6):
            for column_block in range(6):
                cnf.extend(CardEnc.equals([symbol(kind, row_block, column_block, phase) for phase in range(3)], 1, vpool=pool, encoding=EncType.pairwise))

    # Triangle row/column permutations independently make the first repeated
    # outer triple canonical.
    for position, phase in enumerate((0, 1, 2)):
        cnf.append([symbol("r", 0, position, phase)])
        cnf.append([symbol("q", 0, position, phase)])

    def cell(row: int, column: int, phase: int) -> int:
        if row < 3 and column < 3:
            return constant[phase]
        if row < 3:
            column_block = (column - 3) // 2
            return symbol("q", column_block, row, phase)
        if column < 3:
            row_block = (row - 3) // 2
            return symbol("r", row_block, column, phase)
        row_block, row_side = divmod(row - 3, 2)
        column_block, column_side = divmod(column - 3, 2)
        kind = "x" if row_side == column_side else "y"
        return symbol(kind, row_block, column_block, phase)

    def constrain_pair(axis: str, first: int, second: int) -> None:
        differences = [[], [], []]
        for position in range(ORDER):
            for difference in range(3):
                indicator = pool.id((axis, first, second, position, difference))
                differences[difference].append(indicator)
                for phase in range(3):
                    if axis == "row_difference":
                        left = cell(first, position, phase)
                        right = cell(second, position, (phase - difference) % 3)
                    else:
                        left = cell(position, first, phase)
                        right = cell(position, second, (phase - difference) % 3)
                    cnf.append([-left, -right, indicator])
                    cnf.append([-indicator, -left, right])
        for literals, count in zip(differences, target_counts(first, second)):
            cnf.extend(CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter))

    for first in range(ORDER):
        for second in range(first + 1, ORDER):
            constrain_pair("row_difference", first, second)
            constrain_pair("column_difference", first, second)
    return cnf, pool, cell


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", default="glucose42")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    args = parser.parse_args()
    cnf, pool, cell = build()
    with Solver(name=args.solver, bootstrap_with=cnf) as solver:
        timer = Timer(args.timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "model": "reduced Q=81 negative-triangle K3+6K2 block form",
            "solver": args.solver,
            "timeout_seconds": args.timeout_seconds,
            "variables": pool.top,
            "clauses": len(cnf.clauses),
            "result": "sat" if result is True else "unsat" if result is False else "unknown_timeout",
        }
        if result is True:
            model = {literal for literal in solver.get_model() if literal > 0}
            rows = [
                [next(phase for phase in range(3) if cell(row, column, phase) in model) for column in range(ORDER)]
                for row in range(ORDER)
            ]
            payload["rows"] = rows
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
