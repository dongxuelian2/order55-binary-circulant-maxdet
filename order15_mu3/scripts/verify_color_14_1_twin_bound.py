"""Exact twin-leaf PSD obstruction for (14,1,0) energies Q=60 and 69."""

from __future__ import annotations

import json


def certificates() -> list[dict[str, object]]:
    cases = [
        {"total_energy": 60, "internal_energy": 18, "cross_energy": 42, "minimum_row_twins": 10, "minimum_column_active": 5},
        {"total_energy": 69, "internal_energy": 18, "cross_energy": 51, "minimum_row_twins": 9, "minimum_column_active": 4},
        {"total_energy": 69, "internal_energy": 27, "cross_energy": 42, "minimum_row_twins": 10, "minimum_column_active": 4},
    ]
    for case in cases:
        eigenvalue = 15 - case["minimum_row_twins"] * case["minimum_column_active"]
        assert eigenvalue < 0
        case["forced_gram_eigenvalue_upper"] = eigenvalue
    return cases


def main() -> None:
    rows = certificates()
    print(json.dumps({
        "color_partition": [14, 1, 0],
        "cases": rows,
        "mechanism": "Twin leaf differences lie in the opposite Gram-support kernel; the complementary Gram is 15 I_m - v J_m.",
        "theorem": "Total energies Q=60 and Q=69 are impossible.",
    }, indent=2))


if __name__ == "__main__":
    main()
