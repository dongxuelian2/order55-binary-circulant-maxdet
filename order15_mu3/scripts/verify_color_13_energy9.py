"""Exact two-hub Schur enumeration for a 13-row color class of energy 9."""

from __future__ import annotations

import json
from collections import defaultdict
from fractions import Fraction

from maxdet.mu3 import Eisenstein, inner_product_values

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm, valuation
    from order15_mu3.scripts.verify_color_14_1_low_energy import _pair_local_types, reachable_sums
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_bounds import is_rational_eisenstein_norm, valuation
    from verify_color_14_1_low_energy import _pair_local_types, reachable_sums
    from verify_trace_stability import BENCHMARK


def arithmetic_admissible(value: int) -> bool:
    return valuation(value, 3) >= 14 and is_rational_eisenstein_norm(value)


def minimum_cross_losses(catalogue) -> dict[int, Fraction]:
    cross_values = tuple(entry["value"] for entry in catalogue if 0 < entry["norm"] <= 81 and entry["norm"] % 9 == 3)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)
    pair_types = _pair_local_types(cross_values, norm9)
    cross_norms = tuple(sorted({value.norm() for value in cross_values}))
    singles = reachable_sums(cross_norms, 11, 168)
    result = {}
    for h in range(10):
        energy = 39 + 9 * h
        losses = [
            Fraction(numerator, 216) + Fraction(energy - pair_energy, 15)
            for pair_energy, numerator in pair_types
            if energy - pair_energy in singles
        ]
        if losses:
            result[h] = min(losses)
    return result


def main() -> None:
    catalogue = inner_product_values(15)
    norm9 = tuple(entry["value"] for entry in catalogue if entry["norm"] == 9)

    def delta(entry) -> int:
        _a, b, c = entry["counts"][0]
        return (b + 2 * c) % 3

    norm3_by_delta = {
        color_delta: tuple(entry["value"] for entry in catalogue if entry["norm"] == 3 and delta(entry) == color_delta)
        for color_delta in (1, 2)
    }
    assert all(len(values) == 3 for values in norm3_by_delta.values())
    losses = minimum_cross_losses(catalogue)
    leaf_determinant = 216 * 15**11
    assert losses[0] == Fraction(38, 15)
    assert leaf_determinant * (15 - losses[0]) * (15 - losses[1]) < BENCHMARK
    assert all(
        leaf_determinant * (15 - losses[h1]) * (15 - losses[h2]) < BENCHMARK
        for h1 in losses for h2 in losses if h1 + h2 > 0
    )

    candidates = {}
    state_counts = {}
    # Fix the large class to color 0.  In (13,2,0) the two outside rows both
    # have color 1.  In (13,1,1) they have colors 1 and 2.  The orientation of
    # G_{hub1,hub2} then has delta 0 and 2 respectively.
    color_cases = {
        "13_2_0": (1, 1, 0),
        "13_1_1": (1, 2, 2),
    }
    for color_type, (delta_u, delta_v, outside_delta) in color_cases.items():
        norm3_u = norm3_by_delta[delta_u]
        norm3_v = norm3_by_delta[delta_v]
        pair_types: set[tuple[int, int, Eisenstein]] = set()
        for g in (Eisenstein(-3), Eisenstein(3)):
            for u1 in norm3_u:
                for u2 in norm3_u:
                    product_u = u1 * g * u2.conjugate()
                    n_u = 15 * (u1.norm() + u2.norm()) - (2 * product_u.a - product_u.b)
                    for v1 in norm3_v:
                        for v2 in norm3_v:
                            product_v = v1 * g * v2.conjugate()
                            n_v = 15 * (v1.norm() + v2.norm()) - (2 * product_v.a - product_v.b)
                            off_numerator = (
                                15 * u1 * v1.conjugate()
                                - u1 * g * v2.conjugate()
                                - u2 * g.conjugate() * v1.conjugate()
                                + 15 * u2 * v2.conjugate()
                            )
                            pair_types.add((n_u, n_v, off_numerator))

        product_values = {u * v.conjugate() for u in norm3_u for v in norm3_v}
        isolate_sums = {Eisenstein()}
        for _ in range(11):
            isolate_sums = {old + value for old in isolate_sums for value in product_values}

        # Scale every Schur entry by 1080=lcm(216,15).
        corrections: dict[tuple[int, int], set[Eisenstein]] = defaultdict(set)
        for n_u, n_v, off_numerator in pair_types:
            for isolate_sum in isolate_sums:
                corrections[(n_u, n_v)].add(5 * off_numerator + 72 * isolate_sum)

        gram_values = tuple(
            entry["value"] for entry in catalogue
            if entry["norm"] <= 81 and delta(entry) == outside_delta
        )
        determinants = set()
        for (n_u, n_v), scaled_corrections in corrections.items():
            diagonal_u = Fraction(64, 5) - Fraction(n_u, 216)
            diagonal_v = Fraction(64, 5) - Fraction(n_v, 216)
            if leaf_determinant * diagonal_u * diagonal_v <= BENCHMARK:
                continue
            for correction in scaled_corrections:
                for gram_value in gram_values:
                    scaled_off = 1080 * gram_value - correction
                    value = Fraction(leaf_determinant) * (
                        diagonal_u * diagonal_v - Fraction(scaled_off.norm(), 1080**2)
                    )
                    if value > BENCHMARK:
                        assert value.denominator == 1
                        determinants.add(value.numerator)
        candidates[color_type] = {
            "above_record_determinants": len(determinants),
            "arithmetic_survivors": sorted(value for value in determinants if arithmetic_admissible(value)),
        }
        state_counts[color_type] = {
            "pair_local_types": len(pair_types),
            "isolate_product_values": len(product_values),
            "eleven_isolate_product_sums": len(isolate_sums),
            "scaled_correction_states": sum(len(values) for values in corrections.values()),
        }

    print(json.dumps({
        "large_color_class_size": 13,
        "large_class_internal_energy": 9,
        "minimum_cross_loss_by_excess_h": {str(h): str(loss) for h, loss in losses.items()},
        "cross_edge_theorem": "All 26 large-to-small Gram norms must equal 3.",
        "oriented_phase_state_counts": state_counts,
        "color_partition_candidates": candidates,
    }, indent=2))


if __name__ == "__main__":
    main()
