"""Exact certificate that no F7-support Gram matrix beats the record."""

from __future__ import annotations

import json
from fractions import Fraction

from maxdet.mu3 import inner_product_values


BENCHMARK = 2**22 * 3**20 * 19
EXPECTED_FRONTIER = (
    (15, 216, 72),
    (33, 198, 63),
    (42, 189, 54),
    (69, 162, 45),
    (87, 144, 36),
)


def compositions(total: int, parts: int, prefix: tuple[int, ...] = ()):
    if parts == 1:
        yield prefix + (total,)
        return
    for value in range(total + 1):
        yield from compositions(total - value, parts - 1, prefix + (value,))


def main() -> None:
    entries = inner_product_values(15)
    # For quotient counts (a,b,c), delta=sum(exponent differences)=b+2c.
    # The exact congruence q/3=delta^2 mod 3 determines equal/different colors.
    for entry in entries:
        q = int(entry["norm"])
        a, b, c = entry["counts"][0]
        delta = (b + 2 * c) % 3
        assert q % 3 == 0
        assert (q // 3) % 3 == (delta * delta) % 3

    different_color = [entry["value"] for entry in entries if 0 < entry["norm"] <= 81 and entry["norm"] % 9 == 3]
    same_color = [entry["value"] for entry in entries if 0 < entry["norm"] <= 81 and entry["norm"] % 9 == 0]
    local_types: set[tuple[int, int, int]] = set()
    for u in different_color:
        for v in different_color:
            for g in same_color:
                d = 225 - g.norm()
                product = u * g * v.conjugate()
                numerator = 15 * (u.norm() + v.norm()) - (2 * product.a - product.b)
                energy = u.norm() + v.norm() + g.norm()
                local_types.add((energy, d, numerator))

    frontier = []
    for energy, d, numerator in local_types:
        contribution = Fraction(numerator, d)
        dominated = any(
            other_energy <= energy
            and other_d >= d
            and Fraction(other_numerator, other_d) <= contribution
            and (other_energy < energy or other_d > d or Fraction(other_numerator, other_d) < contribution)
            for other_energy, other_d, other_numerator in local_types
        )
        if not dominated:
            frontier.append((energy, d, numerator))
    frontier = sorted(frontier)
    assert tuple(frontier) == EXPECTED_FRONTIER

    profiles = []
    for counts in compositions(7, len(frontier)):
        energy = sum(count * local[0] for count, local in zip(counts, frontier))
        determinant_factor = 1
        schur_loss = Fraction(0)
        for count, (_local_energy, d, numerator) in zip(counts, frontier):
            determinant_factor *= d**count
            schur_loss += count * Fraction(numerator, d)
        determinant_value = Fraction(determinant_factor) * (15 - schur_loss)
        if determinant_value > 0:
            assert determinant_value.denominator == 1
            profiles.append((determinant_value.numerator, energy, counts))
    best = max(profiles)
    assert best[0] == BENCHMARK
    assert best[2] == (7, 0, 0, 0, 0)
    print(json.dumps({
        "local_types_enumerated": len(local_types),
        "pareto_frontier_energy_d_numerator": frontier,
        "seven_block_profiles_checked": len(profiles),
        "maximum_determinant": best[0],
        "maximum_energy": best[1],
        "maximum_profile_counts": best[2],
        "theorem": "A counterexample with F7 nonorthogonality support does not exist.",
    }, indent=2))


if __name__ == "__main__":
    main()
