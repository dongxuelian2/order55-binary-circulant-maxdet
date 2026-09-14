"""Boolean/PB encoding of the all-norm-3 F7 Gram realizability problem."""

from __future__ import annotations

import argparse
import json
from itertools import permutations

from z3 import And, Bool, Or, PbEq, Solver, sat

from maxdet.mu3 import determinant, exponent_matrix, gram_from_exponents


ORDER = 15
PAIRS = tuple((1 + 2 * k, 2 + 2 * k) for k in range(7))


def count_equal(formulas: list[object], value: int) -> object:
    return PbEq([(formula, 1) for formula in formulas], value)


def pair_count_formulas(bits: list[list[list[object]]], first: int, second: int) -> tuple[list[object], ...]:
    result: list[list[object]] = [[], [], []]
    for column in range(ORDER):
        for residue in range(3):
            result[residue].append(Or([
                And(bits[first][column][phase], bits[second][column][(phase - residue) % 3])
                for phase in range(3)
            ]))
    return tuple(result)


def add_exact_counts(solver: Solver, formulas: tuple[list[object], ...], counts: tuple[int, int, int]) -> None:
    for residue in range(3):
        solver.add(count_equal(formulas[residue], counts[residue]))


def norm_three(formulas: tuple[list[object], ...]) -> object:
    return Or([
        And([count_equal(formulas[index], values[index]) for index in range(3)])
        for values in permutations((4, 5, 6))
    ])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-ms", type=int, default=300_000)
    args = parser.parse_args()
    solver = Solver()
    solver.set(timeout=args.timeout_ms)
    bits = [[[Bool(f"b_{row}_{column}_{phase}") for phase in range(3)] for column in range(ORDER)] for row in range(ORDER)]
    for row in range(ORDER):
        for column in range(ORDER):
            solver.add(PbEq([(bits[row][column][phase], 1) for phase in range(3)], 1))

    # Hub column-dephasing and the same row-phase/conjugation/column-order
    # canonicalization used by the integer model.
    canonical_first_leaf = (0,) * 4 + (1,) * 6 + (2,) * 5
    for column in range(ORDER):
        solver.add(bits[0][column][0])
        solver.add(bits[1][column][canonical_first_leaf[column]])

    # Orthogonal leaves have equal exponent sums modulo 3. The cross-pair leaf
    # graph is connected, so all leaves share one sign orbit; conjugation fixes
    # it to the canonical G_0i count triple (4,5,6).
    for leaf in range(2, ORDER):
        formulas = pair_count_formulas(bits, 0, leaf)
        solver.add(And([count_equal(formulas[k], (4, 5, 6)[k]) for k in range(3)]))

    matched = {frozenset(pair) for pair in PAIRS}
    for first in range(1, ORDER):
        for second in range(first + 1, ORDER):
            formulas = pair_count_formulas(bits, first, second)
            if frozenset((first, second)) in matched:
                solver.add(norm_three(formulas))
            else:
                add_exact_counts(solver, formulas, (5, 5, 5))

    result = solver.check()
    payload: dict[str, object] = {
        "encoding": "one-hot Boolean with pseudo-Boolean exact counts",
        "model": "F7 support, all 21 support-edge norms equal 3",
        "timeout_ms": args.timeout_ms,
        "result": str(result),
        "reason_unknown": solver.reason_unknown() if result != sat else "",
    }
    if result == sat:
        model = solver.model()
        exponents = [
            [next(phase for phase in range(3) if model.evaluate(bits[row][column][phase])) for column in range(ORDER)]
            for row in range(ORDER)
        ]
        value = determinant(exponent_matrix(exponents))
        gram = gram_from_exponents(exponents)
        payload.update({
            "determinant": [value.a, value.b],
            "determinant_norm": value.norm(),
            "gram_norms": sorted(gram[i][j].norm() for i in range(ORDER) for j in range(i) if gram[i][j].norm()),
            "matrix": exponents,
        })
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
