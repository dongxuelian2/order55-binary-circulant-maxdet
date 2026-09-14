"""Exact arithmetic and matrix primitives for third roots of unity.

Elements of ``Z[omega]`` are stored as ``a + b*omega`` with
``omega**2 + omega + 1 = 0``.  All certificate-facing operations use Python
integers only.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence


@dataclass(frozen=True, slots=True)
class Eisenstein:
    """An Eisenstein integer ``a + b*omega``."""

    a: int = 0
    b: int = 0

    def __add__(self, other: object) -> "Eisenstein":
        if isinstance(other, int):
            other = Eisenstein(other)
        if not isinstance(other, Eisenstein):
            return NotImplemented
        return Eisenstein(self.a + other.a, self.b + other.b)

    __radd__ = __add__

    def __neg__(self) -> "Eisenstein":
        return Eisenstein(-self.a, -self.b)

    def __sub__(self, other: object) -> "Eisenstein":
        return self + (-_coerce(other))

    def __rsub__(self, other: object) -> "Eisenstein":
        return _coerce(other) - self

    def __mul__(self, other: object) -> "Eisenstein":
        if isinstance(other, int):
            other = Eisenstein(other)
        if not isinstance(other, Eisenstein):
            return NotImplemented
        # omega^2 = -1 - omega
        return Eisenstein(
            self.a * other.a - self.b * other.b,
            self.a * other.b + self.b * other.a - self.b * other.b,
        )

    __rmul__ = __mul__

    def conjugate(self) -> "Eisenstein":
        """Return the complex conjugate, using conjugate(omega)=omega^2."""

        return Eisenstein(self.a - self.b, -self.b)

    def norm(self) -> int:
        """Return ``z * conjugate(z) = a^2 - ab + b^2``."""

        return self.a * self.a - self.a * self.b + self.b * self.b

    def divexact(self, divisor: "Eisenstein") -> "Eisenstein":
        """Divide in ``Z[omega]``, raising when the quotient is not integral."""

        divisor = _coerce(divisor)
        denominator = divisor.norm()
        if denominator == 0:
            raise ZeroDivisionError("division by zero Eisenstein integer")
        numerator = self * divisor.conjugate()
        qa, ra = divmod(numerator.a, denominator)
        qb, rb = divmod(numerator.b, denominator)
        if ra or rb:
            raise ArithmeticError(f"non-exact Eisenstein division: {self}/{divisor}")
        quotient = Eisenstein(qa, qb)
        if quotient * divisor != self:
            raise AssertionError("internal Eisenstein exact-division failure")
        return quotient

    def __complex__(self) -> complex:
        return complex(self.a - self.b / 2, self.b * (3**0.5) / 2)


ZERO = Eisenstein()
ONE = Eisenstein(1)
OMEGA = Eisenstein(0, 1)
OMEGA2 = Eisenstein(-1, -1)
ROOTS = (ONE, OMEGA, OMEGA2)
UNITS = (ONE, -ONE, OMEGA, -OMEGA, OMEGA2, -OMEGA2)


def _coerce(value: object) -> Eisenstein:
    if isinstance(value, Eisenstein):
        return value
    if isinstance(value, int):
        return Eisenstein(value)
    raise TypeError(f"expected Eisenstein integer or int, got {type(value).__name__}")


def exponent_matrix(exponents: Sequence[Sequence[int]]) -> list[list[Eisenstein]]:
    """Convert a rectangular exponent array modulo three to roots of unity."""

    rows = [list(row) for row in exponents]
    if not rows or not rows[0]:
        raise ValueError("matrix must be nonempty")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows must have equal length")
    return [[ROOTS[int(value) % 3] for value in row] for row in rows]


def determinant(matrix: Sequence[Sequence[Eisenstein | int]]) -> Eisenstein:
    """Compute a square determinant by fraction-free Bareiss elimination."""

    work = [[_coerce(value) for value in row] for row in matrix]
    n = len(work)
    if any(len(row) != n for row in work):
        raise ValueError("expected a square matrix")
    if n == 0:
        return ONE
    sign = 1
    previous = ONE
    for k in range(n - 1):
        if work[k][k] == ZERO:
            pivot_row = next((i for i in range(k + 1, n) if work[i][k] != ZERO), None)
            if pivot_row is None:
                return ZERO
            work[k], work[pivot_row] = work[pivot_row], work[k]
            sign = -sign
        pivot = work[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = pivot * work[i][j] - work[i][k] * work[k][j]
                work[i][j] = numerator.divexact(previous)
            work[i][k] = ZERO
        previous = pivot
    return sign * work[-1][-1]


def gram_from_exponents(exponents: Sequence[Sequence[int]]) -> list[list[Eisenstein]]:
    """Return the exact row Gram matrix ``H H*`` of a ternary exponent matrix."""

    rows = [list(row) for row in exponents]
    if not rows or not rows[0]:
        raise ValueError("matrix must be nonempty")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows must have equal length")
    result: list[list[Eisenstein]] = [[ZERO for _ in rows] for _ in rows]
    for i, row_i in enumerate(rows):
        for j in range(i + 1):
            value = sum((ROOTS[(int(x) - int(y)) % 3] for x, y in zip(row_i, rows[j])), ZERO)
            result[i][j] = value
            result[j][i] = value.conjugate()
    return result


def dephase(exponents: Sequence[Sequence[int]]) -> list[list[int]]:
    """Use row/column phase multiplications to zero the first row and column."""

    rows = [list(map(lambda value: int(value) % 3, row)) for row in exponents]
    if not rows or not rows[0]:
        raise ValueError("matrix must be nonempty")
    width = len(rows[0])
    if any(len(row) != width for row in rows):
        raise ValueError("matrix rows must have equal length")
    top = rows[0][:]
    column_adjusted = [[(value - top[j]) % 3 for j, value in enumerate(row)] for row in rows]
    return [[(value - row[0]) % 3 for value in row] for row in column_adjusted]


def inner_product_values(length: int) -> list[dict[str, object]]:
    """Enumerate all sums of ``length`` third roots, grouped by exact value."""

    if length < 0:
        raise ValueError("length must be nonnegative")
    values: dict[Eisenstein, list[tuple[int, int, int]]] = {}
    for count_0 in range(length + 1):
        for count_1 in range(length - count_0 + 1):
            count_2 = length - count_0 - count_1
            value = count_0 * ONE + count_1 * OMEGA + count_2 * OMEGA2
            values.setdefault(value, []).append((count_0, count_1, count_2))
    return [
        {"value": value, "norm": value.norm(), "counts": tuple(counts)}
        for value, counts in sorted(values.items(), key=lambda item: (item[0].norm(), item[0].a, item[0].b))
    ]


def validate_exponents(exponents: Iterable[Iterable[int]], order: int | None = None) -> list[list[int]]:
    """Validate a square matrix with literal entries in ``{0,1,2}``."""

    rows = [list(row) for row in exponents]
    if not rows or any(len(row) != len(rows) for row in rows):
        raise ValueError("expected a nonempty square matrix")
    if order is not None and len(rows) != order:
        raise ValueError(f"expected order {order}, got {len(rows)}")
    if any(type(value) is not int or value not in (0, 1, 2) for row in rows for value in row):
        raise ValueError("exponents must be literal integers in {0,1,2}")
    return rows
