"""Exact simultaneous-Gram SAT test for a canonical top Q=90 phase class."""

from __future__ import annotations

import argparse
import json
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import Eisenstein, determinant, gram_from_exponents


ORDER = 15
TARGET = 289967495625000000
U = Eisenstein(-1, -2)
V_EDGE = Eisenstein(-2, -1)
PRODUCTS = (Eisenstein(-3, 0), Eisenstein(0, -3), Eisenstein(3, 3))
ISOLATE_PRODUCT_COUNTS = (0, 10, 1)
H = Eisenstein(1, -1)


def _v_for_product(product: Eisenstein) -> Eisenstein:
    candidates = (
        Eisenstein(-2, -1), Eisenstein(1, -1), Eisenstein(1, 2),
    )
    return next(value for value in candidates if U * value.conjugate() == product)


def target_gram() -> list[list[Eisenstein]]:
    gram = [[Eisenstein(15 if i == j else 0) for j in range(ORDER)] for i in range(ORDER)]
    gram[0][1] = gram[1][0] = Eisenstein(3)
    for vertex in range(13):
        gram[13][vertex] = U
        gram[vertex][13] = U.conjugate()
    gram[14][0] = gram[14][1] = V_EDGE
    gram[0][14] = gram[1][14] = V_EDGE.conjugate()
    products = [
        product for count, product in zip(ISOLATE_PRODUCT_COUNTS, PRODUCTS)
        for _ in range(count)
    ]
    assert len(products) == 11
    for vertex, product in zip(range(2, 13), products):
        value = _v_for_product(product)
        gram[14][vertex] = value
        gram[vertex][14] = value.conjugate()
    gram[13][14] = H
    gram[14][13] = H.conjugate()
    assert determinant(gram) == Eisenstein(TARGET), determinant(gram)
    return gram


def _counts(value: Eisenstein) -> tuple[int, int, int]:
    third = (15 - value.a - value.b) // 3
    result = (third + value.a, third + value.b, third)
    assert min(result) >= 0 and sum(result) == 15
    return result


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


def build(axes: tuple[str, ...] = ("row", "column")):
    target = target_gram()
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
    if axes == ("row",):
        # Column phases dephase row 0.  Column permutations then put row 1,
        # whose Gram entry with row 0 is +3, into its canonical 7/4/4 form.
        for column in range(ORDER):
            cnf.append([x(0, column, 0)])
            phase = 0 if column < 7 else 1 if column < 11 else 2
            cnf.append([x(1, column, phase)])
    for axis in axes:
        for first in range(ORDER):
            for second in range(first + 1, ORDER):
                differences = [[], [], []]
                for position in range(ORDER):
                    for difference in range(3):
                        indicator = pool.id((axis, first, second, position, difference))
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
                for literals, count in zip(differences, _counts(target[first][second])):
                    cnf.extend(CardEnc.equals(literals, count, vpool=pool, encoding=EncType.seqcounter))
    if axes == ("row",):
        # Ten large-class isolates have identical target neighborhoods.
        for first in range(2, 11):
            _axis_lex(cnf, pool, x, "row", first, first + 1)
        # Residual column permutations preserve the three row-1 phase blocks.
        for start, stop in ((0, 7), (7, 11), (11, 15)):
            for first in range(start, stop - 1):
                _axis_lex(cnf, pool, x, "column", first, first + 1)
    return cnf, pool, x


def solve(
    solver_name: str = "cadical195", timeout_seconds: float = 300,
    axes: tuple[str, ...] = ("row", "column"),
) -> dict[str, object]:
    cnf, pool, x = build(axes)
    with Solver(name=solver_name, bootstrap_with=cnf) as solver:
        timer = Timer(timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "model": "canonical simultaneous top-Q90 row and column Grams",
            "solver": solver_name,
            "variables": pool.top,
            "clauses": len(cnf.clauses),
            "axes": axes,
            "result": "sat" if result is True else "unsat" if result is False else "unknown_timeout",
        }
        if result is True:
            model = {literal for literal in solver.get_model() if literal > 0}
            rows = [
                [next(phase for phase in range(3) if x(row, column, phase) in model) for column in range(ORDER)]
                for row in range(ORDER)
            ]
            if "row" in axes:
                assert gram_from_exponents(rows) == target_gram()
            if "column" in axes:
                assert gram_from_exponents([list(column) for column in zip(*rows)]) == target_gram()
            payload["rows"] = rows
        return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    parser.add_argument("--row-only", action="store_true")
    args = parser.parse_args()
    axes = ("row",) if args.row_only else ("row", "column")
    print(json.dumps(solve(args.solver, args.timeout_seconds, axes), indent=2))


if __name__ == "__main__":
    main()
