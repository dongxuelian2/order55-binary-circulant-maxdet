from fractions import Fraction
from itertools import product
from math import prod

import sympy as sp

from maxdet.order55 import (
    affine_canonical,
    balanced_squares,
    correlations,
    cyclotomic_norms,
    exact55,
    fourth_bound,
    moment_product_upper,
    ryser_global,
    sector_moments,
    sqrt_interval,
)


def determinant(a):
    n = len(a)
    return int(sp.Matrix(n, n, lambda i, j: a[(j - i) % n]).det())


def test_odd_identities_and_pruning_exhaustive_small_orders():
    for n in (3, 5, 7):
        for a in product((0, 1), repeat=n):
            k = sum(a)
            det = determinant(a)
            assert det >= 0
            assert affine_canonical(a) == affine_canonical(a[1:] + a[:1])
            if not 0 < k <= n // 2:
                continue
            comp = determinant(tuple(1 - v for v in a))
            assert k * comp == (n - k) * det
            c = correlations(a)
            half = c[1 : (n + 1) // 2]
            assert 2 * sum(half) == k * (k - 1)
            assert comp <= ryser_global(k, n)
            assert comp <= fourth_bound(k, sum(v * v for v in half), n)
            assert comp <= fourth_bound(
                k, balanced_squares(k * (k - 1) // 2, n // 2), n
            )


def test_moment_bound_rational_grid():
    for values in product(range(1, 5), repeat=4):
        assert prod(values) <= moment_product_upper(
            sum(values), sum(v * v for v in values), 4
        )
    for x in (Fraction(0), Fraction(2), Fraction(9, 4), Fraction(1, 7)):
        lo, hi = sqrt_interval(x)
        assert lo * lo <= x <= hi * hi


def test_order55_exact_paths_and_sector_moments():
    word = "1111011101000011101110011001011100011010010000101000001"
    a = list(map(int, word))
    det, _ = exact55(word)
    assert det == 87835834345428987650257018578111
    assert det == int(
        sp.Matrix(55, 55, lambda i, j: a[(j - i) % 55]).det(method="domain-ge")
    )
    assert det == sum(a) * prod(cyclotomic_norms(word))

    k = sum(a)
    c = correlations(word)
    moments = sector_moments(word)
    assert sum(moments[m][0] for m in (5, 11, 55)) == Fraction(55 * k - k * k, 2)
    assert sum(moments[m][1] for m in (5, 11, 55)) == Fraction(
        55 * sum(v * v for v in c) - k**4, 2
    )
    assert exact55("0" * 55)[0] == exact55("1" * 55)[0] == 0
