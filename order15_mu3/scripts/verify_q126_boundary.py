"""Exact spectral/Schur certificate placing Q=126 below the Q=132 envelope."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_boundary import generic_size13, one_edge_tight
    from order15_mu3.scripts.verify_q117_boundary import schur_two_by_two
    from order15_mu3.scripts.verify_trace_stability import sqrt_lower, stationary_upper
except ModuleNotFoundError:
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_boundary import generic_size13, one_edge_tight
    from verify_q117_boundary import schur_two_by_two
    from verify_trace_stability import sqrt_lower, stationary_upper


def capped_trace_candidates(
    cap: Fraction, variance_energy: int = 252
) -> list[tuple[str, int, int, Fraction]]:
    """KKT candidates at trace 225, fixed variance, and lambda_max<=cap."""

    order = 15
    trace = Fraction(225)
    trace_square = Fraction(15 * 15**2 + variance_energy)
    candidates = []

    # Interior stationary points have two eigenvalues.  A certified lower
    # approximation to the high eigenvalue proves the m=1 point violates the
    # cap; retaining any other possibly infeasible point is harmless.
    for multiplicity in range(1, order):
        t = sqrt_lower(Fraction(variance_energy, order * multiplicity * (order - multiplicity)))
        high = Fraction(15) + (order - multiplicity) * t
        if high > cap:
            continue
        candidates.append(("interior", 0, multiplicity, stationary_upper(variance_energy, multiplicity)))

    # At a nonzero boundary maximum, c eigenvalues equal the cap and the
    # remaining eigenvalues again take at most two values.  Boundary zero
    # eigenvalues give determinant zero and need not be listed.
    for capped in range(1, order):
        remaining = order - capped
        remaining_trace = trace - capped * cap
        remaining_square = trace_square - capped * cap**2
        variance = remaining_square - remaining_trace**2 / remaining
        if variance < 0:
            continue
        if variance == 0:
            value = cap**capped * (remaining_trace / remaining) ** remaining
            if value > 0:
                candidates.append(("cap", capped, 0, value))
            continue
        for multiplicity in range(1, remaining):
            t = sqrt_lower(Fraction(variance, remaining * multiplicity * (remaining - multiplicity)))
            high = remaining_trace / remaining + (remaining - multiplicity) * t
            low = remaining_trace / remaining - multiplicity * t
            if low <= 0 or high > cap:
                continue
            value = cap**capped * high**multiplicity * low ** (remaining - multiplicity)
            candidates.append(("cap", capped, multiplicity, value))
    return candidates


def certificate() -> dict[str, object]:
    q132 = max(stationary_upper(264, multiplicity) for multiplicity in range(1, 15))

    # All-same color: at most fourteen nonzero edges.  Thus the support has
    # clique number at most five.  Weighted Cauchy--Schwarz and
    # Motzkin--Straus give rho(E)^2<=2*126*(1-1/5)=1008/5.
    cap = Fraction(146, 5)
    assert Fraction(1008, 5) < (cap - 15) ** 2
    spectral = capped_trace_candidates(cap)
    spectral_upper = max(row[3] for row in spectral)
    maximizing = max(spectral, key=lambda row: row[3])
    assert maximizing[:3] == ("cap", 1, 1)
    assert spectral_upper < q132

    # The (13,1,1) allocations use the attainable different-color norms
    # 3,12,21,39 through this energy range.
    allocations = []
    for internal in range(9, 46, 9):
        for outside in (3, 12, 21, 39):
            cross = 126 - internal - outside
            if cross >= 78 and cross % 9 == 6:
                allocations.append((internal, cross, outside))
    assert len(allocations) == 13
    tight = {(9, 78, 39), (18, 87, 21), (27, 78, 21), (36, 78, 12)}

    generic_rows = []
    for internal, cross, outside in allocations:
        if (internal, cross, outside) in tight:
            continue
        if internal <= 36:
            upper = generic_size13(internal, cross)
        else:
            # Five internal edges imply clique number at most three and
            # rho(E)^2<=60<8^2, hence lambda_max<23.
            upper = stationary_upper_dimension(13, 45) * (
                Fraction(15) - Fraction(cross, 2 * 23)
            ) ** 2
        assert upper < q132
        generic_rows.append({
            "internal": internal, "cross": cross, "outside_norm": outside,
            "over_q132": float(upper / q132),
        })

    lambda_p3 = Fraction(193, 10)
    tight_bounds = {
        (9, 78, 39): one_edge_tight(39),
        (18, 87, 21): max(
            schur_two_by_two(3105 * 15**10, lambda_p3, 87),
            schur_two_by_two(216**2 * 15**9, Fraction(18), 87),
        ),
        (27, 78, 21): max(
            schur_two_by_two(198 * 15**11, Fraction(21)),
            schur_two_by_two(216**3 * 15**7, Fraction(18)),
            schur_two_by_two(3105 * 216 * 15**8, lambda_p3),
            schur_two_by_two(44631 * 15**9, Fraction(20)),
            schur_two_by_two(44550 * 15**9, Fraction(21)),
            schur_two_by_two(3024 * 15**10, Fraction(21)),
        ),
        (36, 78, 12): stationary_upper_dimension(13, 36) * (
            Fraction(15) - Fraction(39, 22)
        ) ** 2,
    }
    assert all(value < q132 for value in tight_bounds.values())

    return {
        "q132_envelope_numerator": q132.numerator,
        "q132_envelope_denominator": q132.denominator,
        "all_same_lambda_cap": [cap.numerator, cap.denominator],
        "all_same_kkt_candidates": len(spectral),
        "all_same_maximizing_candidate": list(maximizing[:3]),
        "all_same_over_q132": float(spectral_upper / q132),
        "color_13_1_1_generic": generic_rows,
        "color_13_1_1_tight": [
            {"internal": key[0], "cross": key[1], "outside_norm": key[2], "over_q132": float(value / q132)}
            for key, value in tight_bounds.items()
        ],
        "theorem": "Every Q=126 Gram matrix lies below the Q=132 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
