"""Exact odd-order circulant identities and conservative proof bounds.

Search scores must never enter these functions as proof thresholds.
"""
from fractions import Fraction
from math import gcd, isqrt, prod

PRIMES = (2305843009213696591, 2305843009213697141)


def bits(word):
    a = tuple(int(x) for x in word)
    if not a or any(x not in (0, 1) for x in a):
        raise ValueError('nonempty binary word required')
    return a


def correlations(word):
    a = bits(word)
    n = len(a)
    return tuple(sum(a[t] * a[(t+s) % n] for t in range(n)) for s in range(n))


def fold(word, modulus):
    a = bits(word)
    if len(a) % modulus:
        raise ValueError('modulus must divide order')
    return tuple(sum(a[i::modulus]) for i in range(modulus))


def affine_canonical(word):
    a = bits(word)
    n = len(a)
    return ''.join(map(str, min(tuple(a[(u*i+t) % n] for i in range(n))
                               for u in range(n) if gcd(u, n) == 1
                               for t in range(n))))


def ryser_global(k, n=55):
    if n < 3 or n % 2 == 0 or not 1 <= k <= n//2:
        raise ValueError('odd order and lower paired weight required')
    return (n-k) * Fraction(k*(n-k), n-1)**((n-1)//2)


def sqrt_interval(x, digits=50):
    x = Fraction(x)
    if x < 0:
        raise ValueError('negative radicand')
    scale = 10**digits
    lo = isqrt(x.numerator * scale**2 // x.denominator)
    exact = lo**2 * x.denominator == x.numerator * scale**2
    return Fraction(lo, scale), Fraction(lo if exact else lo+1, scale)


def moment_product_upper(total, square_lower, count):
    """Product bound over nonnegative variables with given sum and min squares.

    Contraction to the mean reduces squares and increases product. At the
    required square sum, any positive stationary point has two values.
    Enumerate their multiplicities with outward rational root intervals.
    Boundary points have zero product. No floating point is used.
    """
    total, square_lower = Fraction(total), Fraction(square_lower)
    if count < 1 or total < 0:
        raise ValueError('invalid moment domain')
    if total == 0 or square_lower > total**2:
        return Fraction(0)
    if count == 1:
        return total
    mean = total/count
    variance = max(Fraction(0), square_lower-total**2/count)
    if variance == 0:
        return mean**count
    best = Fraction(0)
    for m in range(1, count):
        low_var = variance * Fraction(count-m, count*m)
        if low_var >= mean**2:
            continue
        high_var = variance * Fraction(m, count*(count-m))
        low = mean-sqrt_interval(low_var)[0]
        high = mean+sqrt_interval(high_var)[1]
        best = max(best, low**m * high**(count-m))
    return best


def balanced_squares(total, count):
    q, r = divmod(total, count)
    return (count-r)*q*q + r*(q+1)**2


def fourth_bound(k, correlation_half_squares, n=55):
    s1 = Fraction(k*(n-k), 2)
    s2 = Fraction(n*(k*k+2*correlation_half_squares)-k**4, 2)
    return (n-k)*moment_product_upper(s1, s2, (n-1)//2)


def modular_determinant(word, p):
    a = bits(word)
    n = len(a)
    if (p-1) % n:
        raise ValueError('order must divide p-1')
    divisors = [d for d in range(2, n+1) if n % d == 0
                and all(d % e for e in range(2, isqrt(d)+1))]
    for b in range(2, p):
        w = pow(b, (p-1)//n, p)
        if all(pow(w, n//d, p) != 1 for d in divisors):
            break
    return prod(sum(x*pow(w, j*t, p) for t, x in enumerate(a)) % p
                for j in range(n)) % p


def exact55(word):
    a = bits(word)
    if len(a) != 55:
        raise ValueError('order 55 required')
    p, q = PRIMES
    # The odd-circulant Ryser bound applies globally via complement duality.
    upper = max(ryser_global(k) for k in range(1, 28))
    if 2*upper >= p*q:
        raise ArithmeticError('insufficient signed CRT capacity')
    r, s = (modular_determinant(a, prime) for prime in PRIMES)
    value = r + p*((s-r)*pow(p, -1, q) % q)
    if value > p*q//2:
        value -= p*q
    if value < 0 or value > upper:
        raise ArithmeticError('CRT reconstruction outside proven bounds')
    return value, (r, s)


def cyclotomic_norms(word):
    import sympy as sp
    a = bits(word)
    if len(a) != 55:
        raise ValueError('order 55 required')
    x = sp.Symbol('x')
    f = sum(v*x**i for i, v in enumerate(a))
    return tuple(int(sp.resultant(f, sp.cyclotomic_poly(m, x), x))
                 for m in (5, 11, 55))


def sector_moments(word):
    """Exact first/second q power sums in root-order 5,11,55 sectors."""
    a = bits(word)
    if len(a) != 55:
        raise ValueError('order 55 required')
    k = sum(a)
    c = correlations(a)
    def folded_moments(m):
        r = fold(a, m)
        h = [sum(c[t] for t in range(i, 55, m)) for i in range(m)]
        return (Fraction(m*sum(v*v for v in r)-k*k, 2),
                Fraction(m*sum(v*v for v in h)-k**4, 2))
    five, eleven = folded_moments(5), folded_moments(11)
    full = (Fraction(55*k-k*k, 2), Fraction(55*sum(v*v for v in c)-k**4, 2))
    primitive = tuple(full[i]-five[i]-eleven[i] for i in (0, 1))
    return {5: five, 11: eleven, 55: primitive}
