"""Enumerate exact phase descriptors attaining the top Q=90 Schur value."""

from __future__ import annotations

import json
from fractions import Fraction

from maxdet.mu3 import Eisenstein, inner_product_values

try:
    from order15_mu3.scripts.verify_bounds import is_rational_eisenstein_norm
    from order15_mu3.scripts.verify_trace_stability import BENCHMARK
except ModuleNotFoundError:
    from verify_bounds import is_rational_eisenstein_norm
    from verify_trace_stability import BENCHMARK


TARGET = 289967495625000000


def certificate(color_type: str = "13_1_1", target: int = TARGET) -> dict[str, object]:
    catalogue = inner_product_values(15)

    def delta(entry) -> int:
        _a, b, c = entry["counts"][0]
        return (b + 2 * c) % 3

    assert color_type in ("13_1_1", "13_2_0")
    u_values = tuple(entry["value"] for entry in catalogue if entry["norm"] == 3 and delta(entry) == 1)
    if color_type == "13_1_1":
        v_values = tuple(entry["value"] for entry in catalogue if entry["norm"] == 3 and delta(entry) == 2)
        h_values = v_values
    else:
        v_values = u_values
        h_values = (Eisenstein(),)
    product_values = tuple(sorted(
        {u * v.conjugate() for u in u_values for v in v_values}, key=lambda z: (z.a, z.b)
    ))
    assert len(product_values) == 3
    leaf_determinant = 216 * 15**11
    descriptors = []
    arithmetic_values = set()
    for g in (Eisenstein(-3), Eisenstein(3)):
        for u0 in u_values:
            for u1 in u_values:
                product_u = u0 * g * u1.conjugate()
                n_u = 15 * (u0.norm() + u1.norm()) - (2 * product_u.a - product_u.b)
                diagonal_u = Fraction(64, 5) - Fraction(n_u, 216)
                for v0 in v_values:
                    for v1 in v_values:
                        product_v = v0 * g * v1.conjugate()
                        n_v = 15 * (v0.norm() + v1.norm()) - (2 * product_v.a - product_v.b)
                        diagonal_v = Fraction(64, 5) - Fraction(n_v, 216)
                        pair_off = (
                            15 * u0 * v0.conjugate()
                            - u0 * g * v1.conjugate()
                            - u1 * g.conjugate() * v0.conjugate()
                            + 15 * u1 * v1.conjugate()
                        )
                        for count0 in range(12):
                            for count1 in range(12 - count0):
                                counts = (count0, count1, 11 - count0 - count1)
                                isolate_sum = sum(
                                    (count * value for count, value in zip(counts, product_values)), Eisenstein()
                                )
                                correction = 5 * pair_off + 72 * isolate_sum
                                for h in h_values:
                                    scaled_off = 1080 * h - correction
                                    value = Fraction(leaf_determinant) * (
                                        diagonal_u * diagonal_v - Fraction(scaled_off.norm(), 1080**2)
                                    )
                                    if (
                                        value.denominator == 1
                                        and value > BENCHMARK
                                        and is_rational_eisenstein_norm(value.numerator)
                                    ):
                                        arithmetic_values.add(value.numerator)
                                    if value == target:
                                        descriptors.append({
                                            "g": [g.a, g.b], "u0": [u0.a, u0.b], "u1": [u1.a, u1.b],
                                            "v0": [v0.a, v0.b], "v1": [v1.a, v1.b],
                                            "isolate_product_counts": counts,
                                            "isolate_products": [[z.a, z.b] for z in product_values],
                                            "h": [h.a, h.b],
                                            "n_u": n_u, "n_v": n_v,
                                            "scaled_off": [scaled_off.a, scaled_off.b],
                                        })
    assert descriptors
    signatures = sorted({
        (item["g"][0], item["g"][1], item["n_u"], item["n_v"],
         *item["isolate_product_counts"], *item["scaled_off"])
        for item in descriptors
    })
    return {
        "target": target,
        "color_type": color_type,
        "raw_phase_descriptors": len(descriptors),
        "coarse_signatures": len(signatures),
        "phase_descriptors": descriptors,
        "sample_descriptors": descriptors[:24],
        "arithmetic_values_above_record": sorted(arithmetic_values),
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
