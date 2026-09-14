"""Linear exact μ3 intertwining test for the canonical top-Q90 Gram."""

from __future__ import annotations

import argparse
import json

from z3 import If, Int, Solver, Sum, sat, unsat

from maxdet.mu3 import Eisenstein

try:
    from order15_mu3.scripts.analyze_q90_top import TARGET, certificate as descriptor_certificate
    from order15_mu3.scripts.sat_q90_top import ORDER, target_gram
except ModuleNotFoundError:
    from analyze_q90_top import TARGET, certificate as descriptor_certificate
    from sat_q90_top import ORDER, target_gram


UNITS = (Eisenstein(1), Eisenstein(0, 1), Eisenstein(-1, -1))


def _selected_product(variable, coefficient: Eisenstein, left: bool):
    values = [coefficient * unit if left else unit * coefficient for unit in UNITS]
    return (
        If(variable == 0, values[0].a, If(variable == 1, values[1].a, values[2].a)),
        If(variable == 0, values[0].b, If(variable == 1, values[1].b, values[2].b)),
    )


def _eisenstein(pair) -> Eisenstein:
    return Eisenstein(pair[0], pair[1])


def top_grams(
    color_type: str = "13_1_1", target: int = TARGET
) -> tuple[list[list[Eisenstein]], ...]:
    result = []
    seen = set()
    for item in descriptor_certificate(color_type, target)["phase_descriptors"]:
        gram = [[Eisenstein(15 if i == j else 0) for j in range(ORDER)] for i in range(ORDER)]
        g = _eisenstein(item["g"])
        u = _eisenstein(item["u0"])
        v_edge = _eisenstein(item["v0"])
        gram[0][1] = g
        gram[1][0] = g.conjugate()
        for vertex in range(13):
            gram[13][vertex] = u
            gram[vertex][13] = u.conjugate()
        gram[14][0] = gram[14][1] = v_edge
        gram[0][14] = gram[1][14] = v_edge.conjugate()
        products = [
            _eisenstein(product)
            for count, product in zip(item["isolate_product_counts"], item["isolate_products"])
            for _ in range(count)
        ]
        for vertex, product in zip(range(2, 13), products):
            quotient_numerator = u.conjugate() * product
            assert quotient_numerator.a % 3 == 0 and quotient_numerator.b % 3 == 0
            value = Eisenstein(
                quotient_numerator.a // 3, quotient_numerator.b // 3
            ).conjugate()
            assert u * value.conjugate() == product
            gram[14][vertex] = value
            gram[vertex][14] = value.conjugate()
        h = _eisenstein(item["h"])
        gram[13][14] = h
        gram[14][13] = h.conjugate()
        key = tuple((value.a, value.b) for row in gram for value in row)
        if key not in seen:
            seen.add(key)
            result.append(gram)
    return tuple(result)


def _canonical_switch_signature(
    gram: list[list[Eisenstein]], allow_conjugation: bool = True
) -> tuple[tuple[int, int], ...]:
    """Canonicalize μ3 diagonal switching, the support automorphisms, and conjugation."""

    signatures = []
    for swap_edge in (False, True):
        for swap_hubs in (False, True):
            for conjugate in ((False, True) if allow_conjugation else (False,)):
                permutation = list(range(ORDER))
                if swap_edge:
                    permutation[0], permutation[1] = permutation[1], permutation[0]
                if swap_hubs:
                    permutation[13], permutation[14] = permutation[14], permutation[13]
                matrix = [[gram[permutation[i]][permutation[j]] for j in range(ORDER)] for i in range(ORDER)]
                if conjugate:
                    matrix = [[value.conjugate() for value in row] for row in matrix]

                phases = [None] * ORDER
                phases[0] = UNITS[0]
                tree = [(0, 1), (0, 13), (0, 14)] + [(13, vertex) for vertex in range(2, 13)]
                for parent, child in tree:
                    orbit = tuple(matrix[parent][child] * unit for unit in UNITS)
                    representative = min(orbit, key=lambda value: (value.a, value.b))
                    phases[child] = next(
                        unit for unit in UNITS
                        if phases[parent] * matrix[parent][child] * unit.conjugate() == representative
                    )
                switched = [
                    [phases[i] * matrix[i][j] * phases[j].conjugate() for j in range(ORDER)]
                    for i in range(ORDER)
                ]
                core = (0, 1, 13, 14)
                signature = [
                    (switched[i][j].a, switched[i][j].b)
                    for i in core for j in core if i < j
                ]
                isolate_signatures = sorted(
                    ((switched[13][i].a, switched[13][i].b),
                     (switched[14][i].a, switched[14][i].b))
                    for i in range(2, 13)
                )
                signature.extend(value for pair in isolate_signatures for value in pair)
                signatures.append(tuple(signature))
    return min(signatures)


def orbit_certificate(color_type: str = "13_1_1", target: int = TARGET) -> dict[str, object]:
    grams = top_grams(color_type, target)
    signatures = {_canonical_switch_signature(gram) for gram in grams}
    switching_signatures = {
        _canonical_switch_signature(gram, allow_conjugation=False) for gram in grams
    }
    assert grams and signatures
    return {
        "color_type": color_type,
        "target": target,
        "top_phase_grams": len(grams),
        "monomial_orbits": len(switching_signatures),
        "monomial_conjugacy_orbits": len(signatures),
        "theorem": f"The top {color_type} phase Grams have the displayed exact orbit counts.",
    }


def orbit_representatives(
    color_type: str = "13_1_1", target: int = TARGET
) -> tuple[list[list[Eisenstein]], ...]:
    representatives = {}
    for gram in top_grams(color_type, target):
        signature = _canonical_switch_signature(gram, allow_conjugation=False)
        representatives.setdefault(signature, gram)
    return tuple(representatives.values())


def intertwining_certificate(
    timeout_seconds: float = 60, color_type: str = "13_1_1", target: int = TARGET
) -> dict[str, object]:
    representatives = orbit_representatives(color_type, target)
    rows = []
    for row_index, row_gram in enumerate(representatives):
        for column_index, column_gram in enumerate(representatives):
            solver, _x = build(False, row_gram, column_gram)
            solver.set(timeout=int(timeout_seconds * 1000))
            result = solver.check()
            rows.append({
                "row_orbit": row_index,
                "column_orbit": column_index,
                "result": "sat" if result == sat else "unsat" if result == unsat else "unknown_timeout",
            })
    assert all(row["result"] == "unsat" for row in rows)
    return {
        **orbit_certificate(color_type, target),
        "pair_results": rows,
        "theorem": f"No top {color_type} row/column Gram pair admits a mu3-valued intertwiner.",
    }
def build(
    include_row_gram: bool = False,
    row_gram: list[list[Eisenstein]] | None = None,
    column_gram: list[list[Eisenstein]] | None = None,
):
    row_gram = target_gram() if row_gram is None else row_gram
    column_gram = row_gram if column_gram is None else column_gram
    solver = Solver()
    x = [[Int(f"x_{row}_{column}") for column in range(ORDER)] for row in range(ORDER)]
    for row in range(ORDER):
        for column in range(ORDER):
            solver.add(x[row][column] >= 0, x[row][column] <= 2)
    solver.add(x[0][0] == 0)

    # G H = H G, two exact Eisenstein-coordinate equations per entry.
    for row in range(ORDER):
        for column in range(ORDER):
            left = [_selected_product(x[k][column], row_gram[row][k], True) for k in range(ORDER)]
            right = [_selected_product(x[row][k], column_gram[k][column], False) for k in range(ORDER)]
            solver.add(Sum(*(item[0] for item in left)) == Sum(*(item[0] for item in right)))
            solver.add(Sum(*(item[1] for item in left)) == Sum(*(item[1] for item in right)))

    if include_row_gram:
        # Exact HH*=G.  A phase difference d selects 1,omega,omega^2.
        for first in range(ORDER):
            for second in range(first + 1, ORDER):
                a_terms = []
                b_terms = []
                for column in range(ORDER):
                    difference = (x[first][column] - x[second][column]) % 3
                    a_terms.append(If(difference == 0, 1, If(difference == 1, 0, -1)))
                    b_terms.append(If(difference == 0, 0, If(difference == 1, 1, -1)))
                solver.add(Sum(*a_terms) == row_gram[first][second].a)
                solver.add(Sum(*b_terms) == row_gram[first][second].b)
    return solver, x


def solve(timeout_seconds: float = 60, include_row_gram: bool = False) -> dict[str, object]:
    solver, x = build(include_row_gram)
    solver.set(timeout=int(timeout_seconds * 1000))
    result = solver.check()
    payload: dict[str, object] = {
        "model": "mu3 matrix commuting with canonical top-Q90 Gram",
        "include_row_gram": include_row_gram,
        "timeout_seconds": timeout_seconds,
        "result": "sat" if result == sat else "unsat" if str(result) == "unsat" else "unknown_timeout",
    }
    if result == sat:
        model = solver.model()
        payload["rows"] = [[model.evaluate(x[i][j]).as_long() for j in range(ORDER)] for i in range(ORDER)]
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout-seconds", type=float, default=60)
    parser.add_argument("--include-row-gram", action="store_true")
    args = parser.parse_args()
    print(json.dumps(solve(args.timeout_seconds, args.include_row_gram), indent=2))


if __name__ == "__main__":
    main()
