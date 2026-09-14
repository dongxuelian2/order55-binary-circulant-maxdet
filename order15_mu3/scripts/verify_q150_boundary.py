"""Exact color/support certificate placing Q=150 below Q=153."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_boundary import generic_size13, one_edge_tight
    from order15_mu3.scripts.verify_q132_boundary import schur_with_residual
    from order15_mu3.scripts.verify_q141_boundary import (
        _energy18_upper, _energy27_upper, _paw_energy36_upper,
        certificate as q141_certificate,
    )
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_boundary import generic_size13, one_edge_tight
    from verify_q132_boundary import schur_with_residual
    from verify_q141_boundary import _energy18_upper, _energy27_upper, _paw_energy36_upper, certificate as q141_certificate
    from verify_trace_stability import stationary_upper


def certificate() -> dict[str, object]:
    q153 = max(stationary_upper(306, multiplicity) for multiplicity in range(1, 15))
    prior = q141_certificate()
    exact = prior["size13_exact_internal_maxima"]
    outside_norms = (0, 9, 27, 36, 63, 81)
    universal = {18: _energy18_upper(), 27: _energy27_upper()}
    caps = {36: Fraction(22), 45: Fraction(23), 54: Fraction(24), 63: Fraction(25), 72: Fraction(127, 5)}
    rows = []
    values = []
    for internal in range(9, 73, 9):
        for outside in outside_norms:
            cross = 150 - internal - outside
            if cross < 78 or cross % 9 != 6:
                continue
            if internal == 9:
                upper = generic_size13(9, cross) if cross >= 105 else one_edge_tight(outside)
                method = "one-edge generic/residual"
            elif internal in universal:
                upper = universal[internal]
                method = "exact local inverse loss"
            elif internal == 36 and cross == 78:
                upper = max(
                    _paw_energy36_upper(outside),
                    schur_with_residual(1674039150000000, Fraction(21), Fraction(9), outside),
                )
                method = "phase-split four-edge residual"
            else:
                block = Fraction(exact[internal]) if internal in exact else stationary_upper_dimension(13, internal)
                upper = block * (Fraction(15) - Fraction(cross, 2) / caps[internal]) ** 2
                method = "internal determinant plus Schur trace"
            assert upper < q153, (internal, cross, outside, float(upper / q153))
            values.append(upper)
            rows.append({
                "internal": internal, "cross": cross, "outside_norm": outside,
                "method": method, "over_q153": float(upper / q153),
            })

    # Q=150 is 6 mod 9, so the other color type is (14,1).  Internally
    # isolated supports are excluded by the equal-sign orbit lemma.  Generic
    # exact trace and weighted clique caps close every no-isolate energy.
    size14_data = (
        (63, 87, Fraction(25)), (72, 78, Fraction(127, 5)),
        (81, 69, Fraction(261, 10)), (90, 60, Fraction(27)),
        (99, 51, Fraction(138, 5)), (108, 42, Fraction(563, 20)),
    )
    size14 = {}
    for internal, cross, cap in size14_data:
        if internal == 108:
            assert (cap - 15) ** 2 > Fraction(864, 5)
        upper = stationary_upper_dimension(14, internal) * (Fraction(15) - Fraction(cross) / cap)
        assert upper < q153
        size14[internal] = upper

    certified_upper = max(*values, *size14.values())
    return {
        "q153_envelope_numerator": q153.numerator,
        "q153_envelope_denominator": q153.denominator,
        "size13_rows": rows,
        "size13_max_over_q153": max(row["over_q153"] for row in rows),
        "size14_over_q153": {str(key): float(value / q153) for key, value in size14.items()},
        "certified_upper_numerator": certified_upper.numerator,
        "certified_upper_denominator": certified_upper.denominator,
        "theorem": "Every Q=150 Gram matrix lies below the Q=153 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
