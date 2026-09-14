"""PySAT search for 13 pairwise orthogonal ternary phase rows of length 15.

Such a set is necessary for a (14,1,0) leaf block with only one nonzero
internal edge: delete either endpoint of that edge.  The encoding fixes the
zero word, row phases, and one balanced word under column permutations.
"""

from __future__ import annotations

import argparse
import json
from collections import defaultdict
from itertools import product
from threading import Timer

from pysat.card import CardEnc, EncType
from pysat.formula import CNF, IDPool
from pysat.solvers import Solver

from maxdet.mu3 import gram_from_exponents


ORDER = 15


def contingency_tables() -> list[tuple[tuple[int, int, int], ...]]:
    """All joint tables for two balanced rows with balanced differences."""

    compositions = [values for values in product(range(6), repeat=3) if sum(values) == 5]
    result = []
    for rows in product(compositions, repeat=3):
        if not all(sum(rows[a][b] for a in range(3)) == 5 for b in range(3)):
            continue
        if not all(sum(rows[a][b] for a in range(3) for b in range(3) if (a - b) % 3 == d) == 5 for d in range(3)):
            continue
        result.append(rows)
    assert len(result) == 21
    return result


def canonical_second_row(table: tuple[tuple[int, int, int], ...]) -> tuple[int, ...]:
    return tuple(phase for first_phase in range(3) for phase in range(3) for _ in range(table[first_phase][phase]))


def contingency_orbits() -> list[list[int]]:
    """Orbits under row phases, conjugation, and permutations of 0,a,b."""

    tables = contingency_tables()
    index = {table: i for i, table in enumerate(tables)}
    # Affine maps on the exponent pair (a,b).  The nontrivial linear maps
    # swap a,b; swap the zero row with a and dephase; and conjugate globally.
    generators = (
        (((0, 1), (1, 0)), (0, 0)),
        (((-1, 0), (-1, 1)), (0, 0)),
        (((-1, 0), (0, -1)), (0, 0)),
        (((1, 0), (0, 1)), (1, 0)),
        (((1, 0), (0, 1)), (0, 1)),
    )

    def transform(table, matrix, shift):
        result = [[0] * 3 for _ in range(3)]
        for a in range(3):
            for b in range(3):
                new_a = (matrix[0][0] * a + matrix[0][1] * b + shift[0]) % 3
                new_b = (matrix[1][0] * a + matrix[1][1] * b + shift[1]) % 3
                result[new_a][new_b] += table[a][b]
        return tuple(tuple(row) for row in result)

    seen = set()
    orbits = []
    for table_index, table in enumerate(tables):
        if table_index in seen:
            continue
        orbit = {table}
        pending = [table]
        while pending:
            current = pending.pop()
            for matrix, shift in generators:
                image = transform(current, matrix, shift)
                assert image in index
                if image not in orbit:
                    orbit.add(image)
                    pending.append(image)
        indices = sorted(index[item] for item in orbit)
        seen.update(indices)
        orbits.append(indices)
    assert sorted(seen) == list(range(21))
    return orbits


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-seconds", type=float, default=300)
    parser.add_argument("--solver", default="cadical195")
    parser.add_argument("--rows", type=int, default=13)
    parser.add_argument("--all-second-tables", action="store_true")
    parser.add_argument("--orbit-representatives", action="store_true")
    parser.add_argument("--lex-rows", action="store_true")
    parser.add_argument("--lex-columns", action="store_true")
    args = parser.parse_args()
    if not 3 <= args.rows <= 15:
        parser.error("--rows must lie between 3 and 15")
    if args.lex_columns and not args.orbit_representatives:
        parser.error("--lex-columns currently requires --orbit-representatives")
    rows_count = args.rows
    all_tables = contingency_tables()
    orbit_representative_indices = [
        next(index for index in orbit if all_tables[index][0][0] > 0)
        for orbit in contingency_orbits()
    ]

    pool = IDPool()
    cnf = CNF()

    def x(row: int, column: int, phase: int) -> int:
        return pool.id(("x", row, column, phase))

    # One-hot cells, with the first row fixed to zero implicitly by omitting
    # it.  Every other row is orthogonal to zero, hence has five occurrences
    # of each phase.  Row phase symmetry fixes its first entry to zero.
    for row in range(1, rows_count):
        for column in range(ORDER):
            cnf.extend(CardEnc.equals([x(row, column, phase) for phase in range(3)], 1, vpool=pool, encoding=EncType.pairwise))
        for phase in range(3):
            cnf.extend(CardEnc.equals([x(row, column, phase) for column in range(ORDER)], 5, vpool=pool, encoding=EncType.seqcounter))
        if not ((args.all_second_tables or args.orbit_representatives) and row == 2):
            cnf.append([x(row, 0, 0)])

    # Column permutations fix the first nonzero row.
    canonical = (0,) * 5 + (1,) * 5 + (2,) * 5
    for column, phase in enumerate(canonical):
        cnf.append([x(1, column, phase)])

    # Pairwise orthogonality means that every phase difference occurs five
    # times.  Difference indicators are linked exactly to the one-hot cells.
    for first in range(1, rows_count):
        for second in range(first + 1, rows_count):
            difference_literals = [[], [], []]
            for column in range(ORDER):
                for difference in range(3):
                    y = pool.id(("d", first, second, column, difference))
                    difference_literals[difference].append(y)
                    for phase in range(3):
                        other_phase = (phase - difference) % 3
                        cnf.append([-x(first, column, phase), -x(second, column, other_phase), y])
                        cnf.append([-y, -x(first, column, phase), x(second, column, other_phase)])
            for difference in range(3):
                cnf.extend(CardEnc.equals(difference_literals[difference], 5, vpool=pool, encoding=EncType.seqcounter))

    if args.lex_rows:
        # Once every row phase has been normalized at column 0, rows 2.. are
        # interchangeable.  Enforce lexicographic order with exact prefix-
        # equality variables.
        first_lex_row = 3 if args.orbit_representatives else 2
        for first in range(first_lex_row, rows_count - 1):
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

    if args.lex_columns:
        # Once rows 0,1,2 are fixed, columns with the same (row1,row2)
        # prefix are interchangeable.  Sort their suffixes.  Each table's
        # clauses are guarded by its branch selector.
        for table_index in orbit_representative_indices:
            table = all_tables[table_index]
            second_row = canonical_second_row(table)
            selector = pool.id(("table_selector", table_index))
            groups = defaultdict(list)
            for column in range(ORDER):
                groups[(canonical[column], second_row[column])].append(column)
            for columns in groups.values():
                for left_column, right_column in zip(columns, columns[1:]):
                    prefix_equal = pool.id(("column_lex_prefix", table_index, left_column, right_column, 3))
                    cnf.append([-selector, prefix_equal])
                    for row in range(3, rows_count):
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

    with Solver(name=args.solver, bootstrap_with=cnf) as solver:
        if args.orbit_representatives:
            # Choose a representative compatible with the row-phase
            # normalization x(2,0)=0.
            table_cases = []
            for representative in orbit_representative_indices:
                table_cases.append((representative, all_tables[representative]))
        elif args.all_second_tables:
            table_cases = list(enumerate(all_tables))
        else:
            table_cases = [(None, None)]
        runs = []
        result = None
        for table_index, table in table_cases:
            assumptions = [] if table is None else [
                x(2, column, phase) for column, phase in enumerate(canonical_second_row(table))
            ]
            if table is not None and args.lex_columns:
                assumptions.append(pool.id(("table_selector", table_index)))
            timer = Timer(args.timeout_seconds, solver.interrupt)
            timer.start()
            result = solver.solve_limited(assumptions=assumptions, expect_interrupt=True)
            timer.cancel()
            runs.append({
                "table_index": table_index if table is not None else None,
                "table": table,
                "result": "sat" if result is True else "unsat" if result is False else "unknown_timeout",
            })
            if result is True:
                break
            if result is None:
                solver.clear_interrupt()
        aggregate_result = True if any(run["result"] == "sat" for run in runs) else (
            False if len(runs) == len(table_cases) and all(run["result"] == "unsat" for run in runs) else None
        )
        result = aggregate_result
        payload: dict[str, object] = {
            "model": f"{rows_count} pairwise orthogonal mu3 rows of length 15",
            "solver": args.solver,
            "timeout_seconds": args.timeout_seconds,
            "variables": pool.top,
            "clauses": len(cnf.clauses),
            "result": "sat" if result is True else "unsat" if result is False else "unknown_timeout",
            "second_table_runs": runs,
            "contingency_orbits": contingency_orbits() if (args.all_second_tables or args.orbit_representatives) else None,
        }
        if result is True:
            true_literals = set(literal for literal in solver.get_model() if literal > 0)
            rows = [[0] * ORDER]
            rows.extend([
                [next(phase for phase in range(3) if x(row, column, phase) in true_literals) for column in range(ORDER)]
                for row in range(1, rows_count)
            ])
            gram = gram_from_exponents(rows)
            assert all(gram[i][j].norm() == 0 for i in range(rows_count) for j in range(i))
            payload["rows"] = rows
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
