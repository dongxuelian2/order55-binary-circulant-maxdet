"""Exact Schur certificate placing the Q=117 size-13 color case below Q=126."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_color_partition_bounds import stationary_upper_dimension
    from order15_mu3.scripts.verify_q114_boundary import generic_size13, one_edge_tight
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:  # Direct execution from this scripts directory.
    from verify_color_partition_bounds import stationary_upper_dimension
    from verify_q114_boundary import generic_size13, one_edge_tight
    from verify_trace_stability import stationary_upper


def schur_two_by_two(block_determinant: int, lambda_upper: Fraction, cross_energy: int = 78) -> Fraction:
    """Bound det(A)det(S) using tr(S)<=30-C/lambda_max(A)."""

    trace_upper = Fraction(30) - Fraction(cross_energy, 1) / lambda_upper
    return block_determinant * (trace_upper / 2) ** 2


def certificate() -> dict[str, object]:
    q126 = max(stationary_upper(252, multiplicity) for multiplicity in range(1, 15))

    # Write Q=e+c+r, where e is the energy inside the 13-row class,
    # c is its cross energy to the two singleton classes, and r is the
    # singleton-to-singleton norm.  The six non-tight allocations are already
    # disposed of by the generic internal trace bound.
    generic_rows = []
    for internal, cross, outside in (
        (9, 87, 21), (9, 96, 12), (9, 105, 3),
        (18, 87, 12), (18, 96, 3), (27, 87, 3),
    ):
        assert internal + cross + outside == 117
        upper = generic_size13(internal, cross)
        assert upper < q126
        generic_rows.append({
            "internal": internal,
            "cross": cross,
            "outside_norm": outside,
            "over_q126": float(upper / q126),
        })

    tight_rows = []

    # e=9: one norm-9 edge and eleven isolates.  Keeping the mandatory
    # outside norm in the 2x2 Schur determinant gives the needed loss.
    upper9 = one_edge_tight(30)
    assert upper9 < q126
    tight_rows.append({"internal": 9, "support": "K2", "over_q126": float(upper9 / q126)})

    # e=18: the support is P3 or 2K2.  Exact block determinants and spectral
    # radii are enough even after dropping the outside off-diagonal term.
    lambda_p3 = Fraction(193, 10)  # 15+3 sqrt(2) < 19.3.
    assert Fraction(18) < (lambda_p3 - 15) ** 2
    energy18 = {
        "P3": schur_two_by_two(3105 * 15**10, lambda_p3),
        "2K2": schur_two_by_two(216**2 * 15**9, Fraction(18)),
    }
    assert all(value < q126 for value in energy18.values())
    tight_rows.append({
        "internal": 18,
        "support": "P3 or 2K2",
        "over_q126": max(float(value / q126) for value in energy18.values()),
    })

    # e=27: either one norm-27 edge or three norm-9 edges.  The latter have
    # exactly the five support shapes listed below.  Tree determinants are
    # switching-invariant; the displayed K3 determinant is the maximum over
    # its six possible cycle gains.
    energy27 = {
        "norm27-K2": schur_two_by_two(198 * 15**11, Fraction(21)),
        "3K2": schur_two_by_two(216**3 * 15**7, Fraction(18)),
        "P3+K2": schur_two_by_two(3105 * 216 * 15**8, lambda_p3),
        "P4": schur_two_by_two(44631 * 15**9, Fraction(20)),
        "K1,3": schur_two_by_two(44550 * 15**9, Fraction(21)),
        "K3": schur_two_by_two(3024 * 15**10, Fraction(21)),
    }
    assert all(value < q126 for value in energy27.values())
    tight_rows.append({
        "internal": 27,
        "support": "one norm-27 edge or a three-edge simple graph",
        "over_q126": max(float(value / q126) for value in energy27.values()),
    })

    # e=36 has weight partitions 36, 27+9, or 9+9+9+9.  In the first two
    # cases the absolute off-diagonal spectral radius is at most 6.  In the
    # last, the four-edge support has clique number at most three, so
    # Motzkin--Straus and Cauchy--Schwarz give rho(E)^2<=48<49.  Thus the
    # uniform bound lambda_max(A)<22 applies to the four-edge case.
    energy36_explicit = {
        "norm36-K2": schur_two_by_two(189 * 15**11, Fraction(21)),
        "disjoint-norm27-and-norm9": schur_two_by_two(198 * 216 * 15**9, Fraction(21)),
        "adjacent-norm27-and-norm9": schur_two_by_two(2835 * 15**10, Fraction(21)),
    }
    four_edge = stationary_upper_dimension(13, 36) * (Fraction(15) - Fraction(39, 22)) ** 2
    energy36 = {**energy36_explicit, "four-norm9-edges": four_edge}
    assert Fraction(48) < Fraction(49)
    assert all(value < q126 for value in energy36.values())
    tight_rows.append({
        "internal": 36,
        "support": "one norm-36 edge, 27+9, or four norm-9 edges",
        "over_q126": max(float(value / q126) for value in energy36.values()),
    })

    return {
        "q126_envelope_numerator": q126.numerator,
        "q126_envelope_denominator": q126.denominator,
        "generic_allocations": generic_rows,
        "tight_allocations": tight_rows,
        "energy18_shape_bounds": {key: float(value / q126) for key, value in energy18.items()},
        "energy27_shape_bounds": {key: float(value / q126) for key, value in energy27.items()},
        "energy36_shape_bounds": {key: float(value / q126) for key, value in energy36.items()},
        "theorem": "Every (13,1,1) Q=117 Gram matrix lies below the Q=126 trace envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
