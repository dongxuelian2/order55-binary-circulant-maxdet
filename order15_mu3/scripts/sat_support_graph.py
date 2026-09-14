"""Exact PySAT realizability search for sparse order-15 mu3 Gram supports.

The first rows are isolated vertices.  Column multiplication fixes row 0 to
zero, column permutations fix row 1, and the five exact contingency-table
orbits fix row 2.  Every nonedge is forced orthogonal; every edge is forced to
have Gram norm 9.  Thus UNSAT for all five branches is an exact obstruction to
the requested support graph (independent of Gram-entry phase enumeration).
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import determinant, gram_from_exponents
try:
    from order15_mu3.scripts.sat_orthogonal_13 import (
        ORDER,
        canonical_second_row,
        contingency_orbits,
        contingency_tables,
    )
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from sat_orthogonal_13 import ORDER, canonical_second_row, contingency_orbits, contingency_tables


SHAPES = {
    # Support vertices follow the isolated prefix and use local labels below.
    "triangle_two_edges": (((0, 1), (0, 2), (1, 2), (3, 4), (5, 6)), 8),
    "cycle5": (((0, 1), (1, 2), (2, 3), (3, 4), (0, 4)), 10),
    "triangle_tail2": (((0, 4), (1, 2), (1, 3), (2, 3), (3, 4)), 10),
}


def add_guarded_equals(cnf: CNF, pool: IDPool, literals: list[int], bound: int, guard: int) -> None:
    encoded = CardEnc.equals(literals, bound, vpool=pool, encoding=EncType.seqcounter)
    cnf.extend([[-guard, *clause] for clause in encoded.clauses])


def add_row_lex(cnf: CNF, pool: IDPool, x, first: int, second: int) -> None:
    prefix_equal = pool.id(("row_lex_prefix", first, second, 0))
    cnf.append([prefix_equal])
    for column in range(ORDER):
        same = pool.id(("row_lex_same", first, second, column))
        for phase in range(3):
            cnf.append([-x(first, column, phase), -x(second, column, phase), same])
            cnf.append([-same, -x(first, column, phase), x(second, column, phase)])
        for left_phase in range(3):
            for right_phase in range(left_phase):
                cnf.append([-prefix_equal, -x(first, column, left_phase), -x(second, column, right_phase)])
        next_prefix = pool.id(("row_lex_prefix", first, second, column + 1))
        cnf.append([-next_prefix, prefix_equal])
        cnf.append([-next_prefix, same])
        cnf.append([-prefix_equal, -same, next_prefix])
        prefix_equal = next_prefix


def add_column_lex(cnf: CNF, pool: IDPool, x, selector: int, table_index: int, second_row) -> None:
    canonical = (0,) * 5 + (1,) * 5 + (2,) * 5
    groups = defaultdict(list)
    for column in range(ORDER):
        groups[(canonical[column], second_row[column])].append(column)
    for columns in groups.values():
        for left_column, right_column in zip(columns, columns[1:]):
            prefix_equal = pool.id(("column_lex_prefix", table_index, left_column, right_column, 3))
            cnf.append([-selector, prefix_equal])
            for row in range(3, ORDER):
                same = pool.id(("column_lex_same", table_index, left_column, right_column, row))
                for phase in range(3):
                    cnf.append([-selector, -x(row, left_column, phase), -x(row, right_column, phase), same])
                    cnf.append([-selector, -same, -x(row, left_column, phase), x(row, right_column, phase)])
                for left_phase in range(3):
                    for right_phase in range(left_phase):
                        cnf.append([-selector, -prefix_equal, -x(row, left_column, left_phase), -x(row, right_column, right_phase)])
                next_prefix = pool.id(("column_lex_prefix", table_index, left_column, right_column, row + 1))
                cnf.append([-selector, -next_prefix, prefix_equal])
                cnf.append([-selector, -next_prefix, same])
                cnf.append([-selector, -prefix_equal, -same, next_prefix])
                prefix_equal = next_prefix


def build(shape_name: str, only_table_index: int | None = None):
    local_edges, isolates = SHAPES[shape_name]
    support_vertices = max(max(edge) for edge in local_edges) + 1
    assert isolates + support_vertices == ORDER
    edges = {tuple(sorted((isolates + a, isolates + b))) for a, b in local_edges}

    pool = IDPool()
    cnf = CNF()

    def x(row: int, column: int, phase: int) -> int:
        return pool.id(("x", row, column, phase))

    # Row 0 is fixed to zero implicitly.  All other rows are orthogonal to it
    # because row 0 is an isolate, so they contain each phase five times.
    for row in range(1, ORDER):
        for column in range(ORDER):
            cnf.extend(CardEnc.equals([x(row, column, phase) for phase in range(3)], 1, vpool=pool, encoding=EncType.pairwise))
        for phase in range(3):
            cnf.extend(CardEnc.equals([x(row, column, phase) for column in range(ORDER)], 5, vpool=pool, encoding=EncType.seqcounter))
        if row != 2:
            cnf.append([x(row, 0, 0)])

    canonical = (0,) * 5 + (1,) * 5 + (2,) * 5
    for column, phase in enumerate(canonical):
        cnf.append([x(1, column, phase)])

    tables = contingency_tables()
    representatives = [next(index for index in orbit if tables[index][0][0] > 0) for orbit in contingency_orbits()]
    if only_table_index is not None:
        if only_table_index not in representatives:
            raise ValueError(f"table index must be one of {representatives}")
        representatives = [only_table_index]
    table_selectors = {index: pool.id(("table_selector", index)) for index in representatives}
    cnf.append(list(table_selectors.values()))
    for first_index, first_selector in table_selectors.items():
        for second_index, second_selector in table_selectors.items():
            if first_index < second_index:
                cnf.append([-first_selector, -second_selector])
    for table_index, selector in table_selectors.items():
        second_row = canonical_second_row(tables[table_index])
        for column, phase in enumerate(second_row):
            cnf.append([-selector, x(2, column, phase)])
        add_column_lex(cnf, pool, x, selector, table_index, second_row)

    # Rows 0,1,2 are distinguished by the contingency-orbit normalization.
    # The remaining isolated rows are interchangeable, but row 2 must not be
    # included in this lex order: the affine map selecting a table orbit need
    # not preserve which remaining isolate is lexicographically first.
    for first in range(3, isolates - 1):
        add_row_lex(cnf, pool, x, first, first + 1)

    norm9_counts = ((7, 4, 4), (4, 7, 4), (4, 4, 7), (6, 6, 3), (6, 3, 6), (3, 6, 6))
    for first in range(1, ORDER):
        for second in range(first + 1, ORDER):
            difference_literals = [[], [], []]
            for column in range(ORDER):
                for difference in range(3):
                    y = pool.id(("d", first, second, column, difference))
                    difference_literals[difference].append(y)
                    for phase in range(3):
                        other_phase = (phase - difference) % 3
                        cnf.append([-x(first, column, phase), -x(second, column, other_phase), y])
                        cnf.append([-y, -x(first, column, phase), x(second, column, other_phase)])
            if (first, second) not in edges:
                for difference in range(3):
                    cnf.extend(CardEnc.equals(difference_literals[difference], 5, vpool=pool, encoding=EncType.seqcounter))
            else:
                selectors = [pool.id(("edge_type", first, second, case)) for case in range(6)]
                cnf.append(selectors)
                for left in range(6):
                    for right in range(left + 1, 6):
                        cnf.append([-selectors[left], -selectors[right]])
                for selector, counts in zip(selectors, norm9_counts):
                    for difference, count in enumerate(counts):
                        add_guarded_equals(cnf, pool, difference_literals[difference], count, selector)
    return cnf, pool, x, edges, representatives, table_selectors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("shape", choices=tuple(SHAPES))
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--timeout-seconds", type=float, default=300)
    parser.add_argument("--table-index", type=int)
    args = parser.parse_args()
    cnf, pool, x, edges, representatives, table_selectors = build(args.shape, args.table_index)
    with Solver(name=args.solver, bootstrap_with=cnf) as solver:
        timer = Timer(args.timeout_seconds, solver.interrupt)
        timer.start()
        result = solver.solve_limited(expect_interrupt=True)
        timer.cancel()
        payload: dict[str, object] = {
            "shape": args.shape,
            "edges": sorted(edges),
            "contingency_representatives": representatives,
            "solver": args.solver,
            "timeout_seconds": args.timeout_seconds,
            "variables": pool.top,
            "clauses": len(cnf.clauses),
            "result": "sat" if result is True else "unsat" if result is False else "unknown_timeout",
        }
        if result is True:
            model = {literal for literal in solver.get_model() if literal > 0}
            rows = [[0] * ORDER] + [
                [next(phase for phase in range(3) if x(row, column, phase) in model) for column in range(ORDER)]
                for row in range(1, ORDER)
            ]
            gram = gram_from_exponents(rows)
            actual_edges = {(i, j) for i in range(ORDER) for j in range(i + 1, ORDER) if gram[i][j].norm() != 0}
            assert actual_edges == edges
            assert all(gram[i][j].norm() == 9 for i, j in edges)
            value = determinant(gram)
            assert value.b == 0
            payload["determinant_norm"] = value.a
            payload["selected_contingency_representative"] = next(index for index, selector in table_selectors.items() if selector in model)
            payload["rows"] = rows
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
