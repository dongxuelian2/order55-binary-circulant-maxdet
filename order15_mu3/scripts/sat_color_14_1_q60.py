"""Exact simultaneous-Gram SAT model for the (14,1,0), Q=60 boundary."""

from __future__ import annotations

import argparse
import json
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import Eisenstein, inner_product_values


ORDER = 15
STAR = Eisenstein(1, -1)
INTERNAL_CASES = {
    "minus": (Eisenstein(-3), Eisenstein(-3)),
    "plus": (Eisenstein(3), Eisenstein(3)),
    "complex_same": (Eisenstein(0, -3), Eisenstein(0, -3)),
    "complex_mixed": (Eisenstein(0, -3), Eisenstein(3, 3)),
    "complex_conjugate": (Eisenstein(3, 3), Eisenstein(3, 3)),
}
DETERMINANT_CLASS = {
    "minus": "minus",
    "plus": "plus",
    "complex_same": "complex",
    "complex_mixed": "complex",
    "complex_conjugate": "complex",
}


def count_lookup() -> dict[Eisenstein, tuple[int, int, int]]:
    return {entry["value"]: entry["counts"][0] for entry in inner_product_values(15)}


def gram_value(first: int, second: int, internal_case: str) -> Eisenstein:
    if second == 14:
        return STAR
    if (first, second) == (0, 1):
        return INTERNAL_CASES[internal_case][0]
    if (first, second) == (2, 3):
        return INTERNAL_CASES[internal_case][1]
    return Eisenstein()


def build(row_case: str, column_case: str):
    if DETERMINANT_CLASS[row_case] != DETERMINANT_CLASS[column_case]:
        raise ValueError("row and column cases must have the same determinant class")
    counts = count_lookup()
    pool = IDPool()
    cnf = CNF()

    def x(row: int, column: int, phase: int) -> int:
        return pool.id(("x", row, column, phase))

    for row in range(ORDER):
        for column in range(ORDER):
            cnf.extend(CardEnc.equals([x(row, column, phase) for phase in range(3)], 1, vpool=pool, encoding=EncType.pairwise))
    cnf.append([x(0, 0, 0)])

    def constrain_pair(axis: str, first: int, second: int, case: str) -> None:
        differences = [[], [], []]
        for position in range(ORDER):
            for difference in range(3):
                indicator = pool.id((axis, first, second, position, difference))
                differences[difference].append(indicator)
                for phase in range(3):
                    if axis == "row_difference":
                        left = x(first, position, phase)
                        right = x(second, position, (phase - difference) % 3)
                    else:
                        left = x(position, first, phase)
                        right = x(position, second, (phase - difference) % 3)
                    cnf.append([-left, -right, indicator])
                    cnf.append([-indicator, -left, right])
        target = counts[gram_value(first, second, case)]
        for literals, count in zip(differences, target):
            cnf.extend(CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter))

    for first in range(ORDER):
        for second in range(first + 1, ORDER):
            constrain_pair("row_difference", first, second, row_case)
            constrain_pair("column_difference", first, second, column_case)

    # Rows 4..13 are the ten singleton leaves and have identical canonical
    # Gram neighborhoods.  Sort them lexicographically; row permutation does
    # not alter the column Gram matrix.
    for first in range(4, 13):
        second = first + 1
        prefix_equal = pool.id(("lex_prefix", first, second, 0))
        cnf.append([prefix_equal])
        for column in range(ORDER):
            same = pool.id(("lex_same", first, second, column))
            for phase in range(3):
                cnf.append([-x(first, column, phase), -x(second, column, phase), same])
                cnf.append([-same, -x(first, column, phase), x(second, column, phase)])
            for left_phase in range(3):
                for right_phase in range(left_phase):
                    cnf.append([-prefix_equal, -x(first, column, left_phase), -x(second, column, right_phase)])
            next_prefix = pool.id(("lex_prefix", first, second, column + 1))
            cnf.append([-next_prefix, prefix_equal])
            cnf.append([-next_prefix, same])
            cnf.append([-prefix_equal, -same, next_prefix])
            prefix_equal = next_prefix
    return cnf, pool, x


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("row_case", choices=tuple(INTERNAL_CASES))
    parser.add_argument("column_case", choices=tuple(INTERNAL_CASES))
    parser.add_argument("--solver", default="glucose42")
    parser.add_argument("--timeout-seconds", type=float, default=180)
    args = parser.parse_args()
    cnf, pool, x = build(args.row_case, args.column_case)
    with Solver(name=args.solver, bootstrap_with=cnf) as solver:
        timer = Timer(args.timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "model": "(14,1,0) Q=60 simultaneous row/column Gram matrices",
            "row_case": args.row_case,
            "column_case": args.column_case,
            "solver": args.solver,
            "variables": pool.top,
            "clauses": len(cnf.clauses),
            "result": "sat" if result is True else "unsat" if result is False else "unknown_timeout",
        }
        if result is True:
            model = {literal for literal in solver.get_model() if literal > 0}
            payload["rows"] = [
                [next(phase for phase in range(3) if x(row, column, phase) in model) for column in range(ORDER)]
                for row in range(ORDER)
            ]
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
