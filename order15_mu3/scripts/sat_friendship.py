"""Exact Z3 realizability test for friendship-graph Gram supports at order 15.

Rows 1..14 are paired. The hub (row 0) is nonorthogonal to every leaf,
paired leaves are nonorthogonal, and leaves in different pairs are orthogonal.
This first model tests the abstract determinant-maximizing profile in which
every support edge has Gram norm 3.
"""

from __future__ import annotations

import argparse
import json
from itertools import permutations

from z3 import If, Int, Or, Solver, Sum, sat

from maxdet.mu3 import determinant, exponent_matrix, gram_from_exponents


ORDER = 15
PAIRS = tuple((1 + 2 * k, 2 + 2 * k) for k in range(7))


def difference_counts(x: list[list[object]], first: int, second: int) -> tuple[object, ...]:
    return tuple(
        Sum([If((x[first][column] - x[second][column]) % 3 == residue, 1, 0) for column in range(ORDER)])
        for residue in range(3)
    )


def norm_three(counts: tuple[object, ...]) -> object:
    return Or([AndList([counts[index] == values[index] for index in range(3)]) for values in permutations((4, 5, 6))])


def AndList(terms: list[object]) -> object:
    # Avoid importing another symbol solely for a three-term conjunction.
    from z3 import And

    return And(terms)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-ms", type=int, default=300_000)
    args = parser.parse_args()
    solver = Solver()
    solver.set(timeout=args.timeout_ms)
    x = [[Int(f"x_{row}_{column}") for column in range(ORDER)] for row in range(ORDER)]
    for row in range(ORDER):
        for column in range(ORDER):
            solver.add(x[row][column] >= 0, x[row][column] <= 2)
    # Column-dephase the hub. Independent leaf row phases rotate each norm-3
    # hub product into one of two sign orbits. Global conjugation fixes the
    # first leaf to the first orbit; residual column permutations sort it.
    for index in range(ORDER):
        solver.add(x[0][index] == 0)
    # G_01 is conjugate to the exponent sum of row 1. The canonical first
    # orbit has exponent counts (4,6,5), and sorting fixes the row completely.
    canonical_first_leaf = (0,) * 4 + (1,) * 6 + (2,) * 5
    for column, value in enumerate(canonical_first_leaf):
        solver.add(x[1][column] == value)

    matched = {frozenset(pair) for pair in PAIRS}
    for leaf in range(2, ORDER):
        counts = difference_counts(x, 0, leaf)
        solver.add(AndList([counts[0] == 4, counts[1] == 5, counts[2] == 6]))
    for first in range(1, ORDER):
        for second in range(first + 1, ORDER):
            counts = difference_counts(x, first, second)
            if frozenset((first, second)) in matched:
                solver.add(norm_three(counts))
            else:
                solver.add(*[counts[residue] == 5 for residue in range(3)])

    result = solver.check()
    payload: dict[str, object] = {
        "model": "F7 support, all 21 support-edge norms equal 3",
        "timeout_ms": args.timeout_ms,
        "result": str(result),
        "reason_unknown": solver.reason_unknown() if result != sat else "",
    }
    if result == sat:
        model = solver.model()
        exponents = [[model.evaluate(x[row][column]).as_long() for column in range(ORDER)] for row in range(ORDER)]
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
