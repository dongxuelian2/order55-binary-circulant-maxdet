"""Exact SAT realizability test for the determinant-maximizing Q=78 paw Gram."""

from __future__ import annotations

import argparse
import json
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import Eisenstein, determinant, gram_from_exponents, inner_product_values


ORDER = 15
HUB = 14
STAR = Eisenstein(1, -1)
PAW_EDGES = ((0, 3), (1, 2), (1, 3), (2, 3))
MAXIMUM = 308697676962890625
LABELINGS = (
    (Eisenstein(0, -3), Eisenstein(0, -3), Eisenstein(3, 0), Eisenstein(3, 3)),
    (Eisenstein(0, -3), Eisenstein(3, 3), Eisenstein(3, 3), Eisenstein(3, 0)),
    (Eisenstein(3, 3), Eisenstein(0, -3), Eisenstein(0, -3), Eisenstein(3, 0)),
    (Eisenstein(3, 3), Eisenstein(3, 3), Eisenstein(3, 0), Eisenstein(0, -3)),
)


def count_lookup() -> dict[Eisenstein, tuple[int, int, int]]:
    return {entry["value"]: entry["counts"][0] for entry in inner_product_values(15)}


def target_gram(first: int, second: int, labeling_index: int) -> Eisenstein:
    if second == HUB:
        return STAR
    edge_values = dict(zip(PAW_EDGES, LABELINGS[labeling_index]))
    return edge_values.get((first, second), Eisenstein())


def add_row_lex(cnf: CNF, pool: IDPool, x, first: int, second: int) -> None:
    prefix = pool.id(("row_lex_prefix", first, second, 0))
    cnf.append([prefix])
    for column in range(ORDER):
        same = pool.id(("row_lex_same", first, second, column))
        for phase in range(3):
            cnf.append([-x(first, column, phase), -x(second, column, phase), same])
            cnf.append([-same, -x(first, column, phase), x(second, column, phase)])
        for left_phase in range(3):
            for right_phase in range(left_phase):
                cnf.append([-prefix, -x(first, column, left_phase), -x(second, column, right_phase)])
        next_prefix = pool.id(("row_lex_prefix", first, second, column + 1))
        cnf.append([-next_prefix, prefix])
        cnf.append([-next_prefix, same])
        cnf.append([-prefix, -same, next_prefix])
        prefix = next_prefix


def add_column_lex(cnf: CNF, pool: IDPool, x, first: int, second: int) -> None:
    """Sort two columns lexicographically; their fixed row-4 phases agree."""

    prefix = pool.id(("column_lex_prefix", first, second, 0))
    cnf.append([prefix])
    for row in range(ORDER):
        same = pool.id(("column_lex_same", first, second, row))
        for phase in range(3):
            cnf.append([-x(row, first, phase), -x(row, second, phase), same])
            cnf.append([-same, -x(row, first, phase), x(row, second, phase)])
        for left_phase in range(3):
            for right_phase in range(left_phase):
                cnf.append([-prefix, -x(row, first, left_phase), -x(row, second, right_phase)])
        next_prefix = pool.id(("column_lex_prefix", first, second, row + 1))
        cnf.append([-next_prefix, prefix])
        cnf.append([-next_prefix, same])
        cnf.append([-prefix, -same, next_prefix])
        prefix = next_prefix


def build(labeling_index: int):
    if not 0 <= labeling_index < len(LABELINGS):
        raise ValueError("labeling index must be 0, 1, 2, or 3")
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

    # Column switching dephases the hub row.  Permuting columns then puts one
    # of the ten interchangeable isolated paw vertices in canonical form.
    for column in range(ORDER):
        cnf.append([x(HUB, column, 0)])
    canonical_isolate = (0,) * 6 + (1,) * 4 + (2,) * 5
    for column, phase in enumerate(canonical_isolate):
        cnf.append([x(4, column, phase)])
    for first in range(4, 13):
        add_row_lex(cnf, pool, x, first, first + 1)
    for group in (range(0, 6), range(6, 10), range(10, 15)):
        for first, second in zip(group, tuple(group)[1:]):
            add_column_lex(cnf, pool, x, first, second)

    for first in range(ORDER):
        for second in range(first + 1, ORDER):
            differences = [[], [], []]
            for column in range(ORDER):
                for difference in range(3):
                    indicator = pool.id(("difference", first, second, column, difference))
                    differences[difference].append(indicator)
                    for phase in range(3):
                        left = x(first, column, phase)
                        right = x(second, column, (phase - difference) % 3)
                        cnf.append([-left, -right, indicator])
                        cnf.append([-indicator, -left, right])
            target = counts[target_gram(first, second, labeling_index)]
            for literals, count in zip(differences, target):
                cnf.extend(CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter))
    return cnf, pool, x


def solve(labeling_index: int, solver_name: str, timeout_seconds: float) -> dict[str, object]:
    cnf, pool, x = build(labeling_index)
    with Solver(name=solver_name, bootstrap_with=cnf) as solver:
        timer = Timer(timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "model": "Q=78 (14,1) maximum paw row Gram",
            "labeling_index": labeling_index,
            "edge_labels": [[value.a, value.b] for value in LABELINGS[labeling_index]],
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
            gram = gram_from_exponents(rows)
            for first in range(ORDER):
                for second in range(first + 1, ORDER):
                    assert gram[first][second] == target_gram(first, second, labeling_index)
            value = determinant(gram)
            assert value == Eisenstein(MAXIMUM)
            payload["rows"] = rows
            payload["determinant"] = value.a
        return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("labeling_index", type=int, choices=range(4))
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    args = parser.parse_args()
    print(json.dumps(solve(args.labeling_index, args.solver, args.timeout_seconds), indent=2))


if __name__ == "__main__":
    main()
