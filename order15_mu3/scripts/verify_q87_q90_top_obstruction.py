"""Exact orbit/intertwining obstruction for the top energy-9 determinant."""

from __future__ import annotations

import json
from functools import lru_cache

try:
    from order15_mu3.scripts.z3_q90_top_commuting import intertwining_certificate
except ModuleNotFoundError:
    from z3_q90_top_commuting import intertwining_certificate


TOP = 289967495625000000
NEXT = 289220156718750000


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    cases = {
        color_type: intertwining_certificate(60, color_type)
        for color_type in ("13_2_0", "13_1_1")
    }
    assert all(
        all(row["result"] == "unsat" for row in result["pair_results"])
        for result in cases.values()
    )
    return {
        "excluded_determinant": TOP,
        "next_arithmetic_energy9_upper": NEXT,
        "color_cases": cases,
        "theorem": "The top energy-9 determinant cannot be realized by an order-15 mu3 matrix.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
