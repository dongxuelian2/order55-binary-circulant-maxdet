"""Exact color-orbit obstruction excluding the entire Q=78 slice."""

from __future__ import annotations

import json

from maxdet.mu3 import Eisenstein, OMEGA, OMEGA2, ONE


ORDER = 15
DELTA = Eisenstein(1, -1)
ROOTS = (ONE, OMEGA, OMEGA2)


def certificate() -> dict[str, object]:
    # At Q=78, congruence excludes the all-same and (13,1,1) partitions.
    # A (13,2) partition costs at least 78 in cross energy and another 9
    # internally in its 13-class, so both Gram color partitions are (14,1).
    assert 78 < 78 + 9

    # Its cross energy is at least 14*3=42.  Hence its internal energy is at
    # most 36 and the 14-vertex internal support has at most four edges.  It
    # therefore has an isolated vertex on each of the row and column sides.
    max_internal_edges = (78 - 14 * 3) // 9
    min_internal_isolates = 14 - 2 * max_internal_edges
    assert max_internal_edges == 4
    assert min_internal_isolates == 6

    # The product of all matrix entries, computed by rows or by columns,
    # shows that the two singleton-minus-large color differences agree.
    # Switching normalizes all cross entries on both sides to the same
    # norm-3 representative delta (or simultaneously to its negative).
    # At an internally isolated row i and column j, (G-15I)H=H(K-15I)
    # reduces to delta*u = v*conj(delta), with u,v third roots.  The two
    # displayed orbits are disjoint, so that equation has no solution.
    left_orbit = {DELTA * root for root in ROOTS}
    right_orbit = {root * DELTA.conjugate() for root in ROOTS}
    assert left_orbit.isdisjoint(right_orbit)
    # Replacing delta by -delta negates both sides and changes nothing.
    assert {-value for value in left_orbit}.isdisjoint({-value for value in right_orbit})

    def pairs(values):
        return sorted([[value.a, value.b] for value in values])

    return {
        "energy": 78,
        "forced_row_and_column_partition": [14, 1, 0],
        "maximum_internal_edges": max_internal_edges,
        "minimum_internal_isolates_per_side": min_internal_isolates,
        "delta_root_orbit": pairs(left_orbit),
        "conjugate_delta_root_orbit": pairs(right_orbit),
        "orbit_intersection": [],
        "theorem": "No order-15 mu3 Gram matrix has total energy Q=78.",
    }


def main() -> None:
    print(json.dumps(certificate(), indent=2))


if __name__ == "__main__":
    main()
