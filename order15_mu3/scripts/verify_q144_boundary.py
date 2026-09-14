"""Exact spectral/Schur certificate placing Q=144 below Q=150."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_boundary import generic_size13, one_edge_tight
    from order15_mu3.scripts.verify_q126_boundary import capped_trace_candidates
    from order15_mu3.scripts.verify_q141_boundary import (
        _energy18_upper,
        _energy27_upper,
        certificate as q141_certificate,
    )
    from order15_mu3.scripts.verify_trace_stability import sqrt_lower, stationary_upper
except ModuleNotFoundError:
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_boundary import generic_size13, one_edge_tight
    from verify_q126_boundary import capped_trace_candidates
    from verify_q141_boundary import _energy18_upper, _energy27_upper, certificate as q141_certificate
    from verify_trace_stability import sqrt_lower, stationary_upper


def one_edge_cross96_outside39() -> Fraction:
    """Schur bound for e=9,c=96,r=39 in the (13,1,1) split."""

    block_determinant = 216 * 15**11
    trace_upper = Fraction(30) - Fraction(96, 18)
    # The two cross-vector energies are 39+57 or 48+48.  Since the one-edge
    # block has lambda_min=12, the inverse cross correction is at most 4.
    correction_upper = Fraction(4)
    residual_lower = sqrt_lower(Fraction(39)) - correction_upper
    assert residual_lower > 0
    return block_determinant * ((trace_upper / 2) ** 2 - residual_lower**2)


def certificate() -> dict[str, object]:
    q150 = max(stationary_upper(300, multiplicity) for multiplicity in range(1, 15))

    # All-same Q=144 has at most sixteen support edges and hence clique number
    # at most six.  Weighted Motzkin--Straus gives rho(E)^2<=240<15.5^2.
    cap = Fraction(61, 2)
    assert Fraction(240) < (cap - 15) ** 2
    all_same = max(row[3] for row in capped_trace_candidates(cap, 288))
    assert all_same < q150

    prior = q141_certificate()
    exact_internal = prior["size13_exact_internal_maxima"]
    universal = {18: _energy18_upper(), 27: _energy27_upper()}
    lambda_caps = {36: Fraction(22), 45: Fraction(23), 54: Fraction(24), 63: Fraction(25)}
    outside_norms = (3, 12, 21, 39, 48, 57)
    rows = []
    values = []
    for internal in range(9, 64, 9):
        for outside in outside_norms:
            cross = 144 - internal - outside
            if cross < 78 or cross % 9 != 6:
                continue
            if internal == 9:
                if cross >= 114:
                    upper = generic_size13(9, cross)
                    method = "generic"
                elif cross == 96:
                    upper = one_edge_cross96_outside39()
                    method = "cross96 residual"
                else:
                    upper = one_edge_tight(outside)
                    method = "minimal-cross residual"
            elif internal in universal:
                upper = universal[internal]
                method = "exact local inverse loss"
            else:
                block = (
                    Fraction(exact_internal[internal])
                    if internal in exact_internal
                    else stationary_upper_dimension(13, internal)
                )
                upper = block * (
                    Fraction(15) - Fraction(cross, 2) / lambda_caps[internal]
                ) ** 2
                method = "internal determinant plus Schur trace"
            assert upper < q150, (internal, cross, outside, method, float(upper / q150))
            values.append(upper)
            rows.append({
                "internal": internal,
                "cross": cross,
                "outside_norm": outside,
                "method": method,
                "over_q150": float(upper / q150),
            })

    certified_upper = max(all_same, *values)
    return {
        "q150_envelope_numerator": q150.numerator,
        "q150_envelope_denominator": q150.denominator,
        "all_same_lambda_cap": [cap.numerator, cap.denominator],
        "all_same_over_q150": float(all_same / q150),
        "size13_rows": rows,
        "size13_max_over_q150": max(row["over_q150"] for row in rows),
        "certified_upper_numerator": certified_upper.numerator,
        "certified_upper_denominator": certified_upper.denominator,
        "theorem": "Every Q=144 Gram matrix lies below the Q=150 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
