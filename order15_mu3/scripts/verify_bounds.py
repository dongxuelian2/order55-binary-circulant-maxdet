"""Replay the finite arithmetic checks in the first order-15 upper bound."""

from __future__ import annotations

import json
from collections import Counter

from sympy import factorint

from maxdet.mu3 import inner_product_values


EXPECTED_SPECTRUM = {
    0: 1, 3: 6, 9: 6, 12: 6, 21: 12, 27: 6, 36: 6,
    39: 12, 48: 6, 57: 12, 63: 12, 75: 6, 81: 3, 84: 6,
    93: 6, 111: 6, 117: 6, 144: 3, 147: 6, 183: 6, 225: 3,
}


def is_rational_eisenstein_norm(value: int) -> bool:
    """Test the norm criterion for a positive rational integer."""

    return value > 0 and all(exponent % 2 == 0 for prime, exponent in factorint(value).items() if prime % 3 == 2)


def valuation(value: int, prime: int) -> int:
    result = 0
    while value % prime == 0:
        value //= prime
        result += 1
    return result


def main() -> None:
    catalogue = inner_product_values(15)
    spectrum = Counter(int(entry["norm"]) for entry in catalogue)
    assert dict(sorted(spectrum.items())) == EXPECTED_SPECTRUM
    edge_norms = sorted(norm for norm in spectrum if 0 < norm < 225)

    # One edge, or two disconnected edges, leaves the odd 5-adic contribution
    # from 15^13 or 15^11 unchanged because every factor 225-q contributes an
    # even valuation.
    assert all(valuation(225 - norm, 5) % 2 == 0 for norm in edge_norms)

    small_sums = sorted({q1 + q2 for q1 in edge_norms for q2 in edge_norms if q1 + q2 < 30})
    assert small_sums == [6, 12, 15, 18, 21, 24]
    path_factors = {total: 3375 - 15 * total for total in small_sums}
    assert all(not is_rational_eisenstein_norm(value) for value in path_factors.values())
    assert is_rational_eisenstein_norm(2925)

    # Triangle: 3240 + 6 sqrt(3) < 3251 follows after squaring 108 < 121.
    assert 108 < 121
    assert not is_rational_eisenstein_norm(3250)

    upper = 3249 * 15**12
    benchmark = 2**22 * 3**20 * 19
    assert 222**3 * 15**9 < upper
    assert benchmark < upper < 15**15
    pair_cutoff = 225 * (1 - benchmark / 15**15)
    assert 81 < pair_cutoff < 84
    print(json.dumps({
        "allowed_norm_spectrum": dict(sorted(spectrum.items())),
        "path_factors_checked": path_factors,
        "rigorous_upper_bound": upper,
        "benchmark": benchmark,
        "squared_gap_factor": upper / benchmark,
        "improvement_pair_norm_cutoff": pair_cutoff,
    }, indent=2))


if __name__ == "__main__":
    main()
