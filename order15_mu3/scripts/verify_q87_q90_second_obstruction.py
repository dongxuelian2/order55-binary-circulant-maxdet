"""Exact orbit/intertwining obstruction for the second energy-9 determinant."""

from __future__ import annotations

import json
from functools import lru_cache

from z3 import sat, unsat

try:
    from order15_mu3.scripts.analyze_q90_top import certificate as descriptor_certificate
    from order15_mu3.scripts.z3_q90_top_commuting import build, orbit_certificate, orbit_representatives
except ModuleNotFoundError:
    from analyze_q90_top import certificate as descriptor_certificate
    from z3_q90_top_commuting import build, orbit_certificate, orbit_representatives


TARGET = 289220156718750000
NEXT = 286464344501953125


def _case_certificate(color_type: str) -> dict[str, object]:
    representatives = orbit_representatives(color_type, TARGET)
    rows = []
    for row_index, row_gram in enumerate(representatives):
        for column_index, column_gram in enumerate(representatives):
            # The sole hard linear pair is the first self-pair in the
            # (13,1,1) case.  Add HH*=G there.  Since det(G)>0 and GH=HK,
            # this also forces H*H=K, so the strengthening is necessary.
            include_row_gram = color_type == "13_1_1" and row_index == column_index == 0
            solver, _x = build(include_row_gram, row_gram, column_gram)
            solver.set(timeout=120000 if include_row_gram else 60000)
            result = solver.check()
            rows.append({
                "row_orbit": row_index,
                "column_orbit": column_index,
                "constraint": "HH*=G and GH=HK" if include_row_gram else "GH=HK",
                "result": "sat" if result == sat else "unsat" if result == unsat else "unknown_timeout",
            })
    assert all(row["result"] == "unsat" for row in rows)
    return {**orbit_certificate(color_type, TARGET), "pair_results": rows}


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    descriptors = {
        color_type: descriptor_certificate(color_type, TARGET)
        for color_type in ("13_2_0", "13_1_1")
    }
    arithmetic_values = sorted({
        value
        for result in descriptors.values()
        for value in result["arithmetic_values_above_record"]
    })
    assert arithmetic_values[-2:] == [TARGET, 289967495625000000]
    assert max(value for value in arithmetic_values if value < TARGET) == NEXT
    cases = {color_type: _case_certificate(color_type) for color_type in descriptors}
    return {
        "excluded_determinant": TARGET,
        "next_arithmetic_energy9_upper": NEXT,
        "color_cases": cases,
        "theorem": "The second energy-9 determinant cannot be realized by an order-15 mu3 matrix.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
