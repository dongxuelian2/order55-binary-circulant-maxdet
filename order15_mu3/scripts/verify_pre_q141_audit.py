"""Exact audit of the all-same Q=126 case against Q=141."""

from __future__ import annotations

import json
from fractions import Fraction

try:
    from order15_mu3.scripts.verify_q135_boundary import fixed_isolate_capped_upper
    from order15_mu3.scripts.verify_trace_stability import stationary_upper
except ModuleNotFoundError:
    from verify_q135_boundary import fixed_isolate_capped_upper
    from verify_trace_stability import stationary_upper


def certificate() -> dict[str, object]:
    q141 = max(stationary_upper(282, multiplicity) for multiplicity in range(1, 15))

    # If the support has no K5, omega<=4 gives rho(E)^2<=189<13.8^2.
    no_k5 = max(row[3] for row in __import__(
        "order15_mu3.scripts.verify_q126_boundary", fromlist=["capped_trace_candidates"]
    ).capped_trace_candidates(Fraction(144, 5), 252))
    assert no_k5 < q141

    # A K5 needs ten edges.  The weighted 10-, 11-, and 12-edge partitions
    # force at least ten, eight, and six isolates; the global 29.2 cap closes
    # all three.  The only 14-edge partition is pure unit weight.  If r of
    # its four extra edges join the K5 component, at least 2+r vertices are
    # isolated.  Motzkin--Straus on that component gives the displayed caps.
    weighted = {
        edges: fixed_isolate_capped_upper(isolates, 252, Fraction(146, 5))
        for edges, isolates in ((10, 10), (11, 8), (12, 6))
    }
    pure_caps = (Fraction(27), Fraction(138, 5), Fraction(141, 5), Fraction(287, 10), Fraction(146, 5))
    pure = {
        r: fixed_isolate_capped_upper(2 + r, 252, cap)
        for r, cap in enumerate(pure_caps)
    }
    assert all(value < q141 for value in weighted.values())
    assert all(value < q141 for value in pure.values())
    return {
        "q141_envelope_numerator": q141.numerator,
        "q141_envelope_denominator": q141.denominator,
        "no_k5_over_q141": float(no_k5 / q141),
        "weighted_k5_over_q141": {str(key): float(value / q141) for key, value in weighted.items()},
        "pure_k5_tradeoff_over_q141": {str(key): float(value / q141) for key, value in pure.items()},
        "theorem": "The all-same Q=126 case lies below the Q=141 envelope.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
