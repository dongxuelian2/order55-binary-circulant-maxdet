"""Exact one-vector Schur losses for sparse 13-row internal blocks."""

from __future__ import annotations

import json
from fractions import Fraction
from functools import lru_cache
from itertools import product

from maxdet.mu3 import Eisenstein, inner_product_values

try:
    from order15_mu3.scripts.verify_color_14_1_energy27 import component_types, orbit_representatives
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_color_14_1_energy27 import component_types, orbit_representatives
    from verify_trace_stability import stationary_upper


EXPECTED_LOCAL_LOSSES = {
    "K2-q9": {216: Fraction(1, 3)},
    "K2-q27": {198: Fraction(7, 22)},
    "P3": {3105: Fraction(11, 23)},
    "P4": {44631: Fraction(18, 29)},
    "K1,3": {44550: Fraction(7, 11)},
    "K3": {
        2916: Fraction(1, 2),
        2943: Fraction(53, 109),
        2997: Fraction(17, 37),
        3024: Fraction(3, 7),
    },
}


def _minimum_losses(states) -> dict[int, Fraction]:
    result = {}
    for _cross_energy, block_determinant, numerator in states:
        loss = Fraction(numerator, block_determinant)
        result[block_determinant] = min(result.get(block_determinant, loss), loss)
    return result


@lru_cache(maxsize=1)
def certificate() -> dict[str, object]:
    catalogue = inner_product_values(15)

    def delta(entry) -> int:
        _a, b, c = entry["counts"][0]
        return (b + 2 * c) % 3

    cross = tuple(entry["value"] for entry in catalogue if entry["norm"] == 3 and delta(entry) == 1)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    norm27 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 27)
    signs = (Eisenstein(-3), Eisenstein(3))
    shapes = {
        "K2-q9": component_types(2, ((0, 1),), (signs,), cross),
        "K2-q27": component_types(2, ((0, 1),), (orbit_representatives(norm27),), cross),
        "P3": component_types(3, ((0, 1), (1, 2)), (signs, signs), cross),
        "P4": component_types(4, ((0, 1), (1, 2), (2, 3)), (signs,) * 3, cross),
        "K1,3": component_types(4, ((0, 1), (0, 2), (0, 3)), (signs,) * 3, cross),
        "K3": component_types(3, ((0, 1), (1, 2), (0, 2)), (signs, signs, norm9), cross),
    }
    losses = {name: _minimum_losses(states) for name, states in shapes.items()}
    assert losses == EXPECTED_LOCAL_LOSSES

    # Each of the two hub vectors has norm-energy 39.  Isolated internal
    # coordinates contribute 3/15=1/5 apiece.  Dropping the Schur
    # off-diagonal gives det(S)<=S11*S22.
    bounds18 = {
        "P3": 3105 * 15**10 * (Fraction(15) - (Fraction(2) + losses["P3"][3105])) ** 2,
        "2K2": 216**2 * 15**9 * (
            Fraction(15) - (Fraction(9, 5) + 2 * losses["K2-q9"][216])
        ) ** 2,
    }

    bounds27 = {}
    support_data = {
        "norm27-K2": (["K2-q27"], 11),
        "3K2": (["K2-q9"] * 3, 7),
        "P3+K2": (["P3", "K2-q9"], 8),
        "P4": (["P4"], 9),
        "K1,3": (["K1,3"], 9),
        "K3": (["K3"], 10),
    }
    for name, (components, isolates) in support_data.items():
        for choices in product(*(losses[component].items() for component in components)):
            determinant_factor = 15**isolates
            loss = Fraction(isolates, 5)
            for block_determinant, local_loss in choices:
                determinant_factor *= block_determinant
                loss += local_loss
            bounds27[f"{name}-{','.join(str(choice[0]) for choice in choices)}"] = (
                determinant_factor * (Fraction(15) - loss) ** 2
            )

    q135 = max(stationary_upper(270, multiplicity) for multiplicity in range(1, 15))
    assert all(value < q135 for value in bounds18.values())
    assert all(value < q135 for value in bounds27.values())
    return {
        "local_minimum_losses": {
            name: {str(det): str(loss) for det, loss in values.items()} for name, values in losses.items()
        },
        "energy18_bounds_over_q135": {name: float(value / q135) for name, value in bounds18.items()},
        "energy27_bounds_over_q135": {name: float(value / q135) for name, value in bounds27.items()},
        "theorem": "Every minimal-cross energy-18 or energy-27 size-13 Schur case lies below Q=135.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
