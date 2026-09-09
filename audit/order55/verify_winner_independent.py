"""Clean-room exact checks for the order-55 winner.

This file deliberately does not import ``maxdet`` or any production script.
It has three separate determinant paths:

* fraction-free Bareiss over Z;
* modular Gaussian elimination at primes not used by production;
* a cyclotomic resultant product using SymPy's public resultant primitive.

The Fourier path uses independently selected primes p == 1 (mod 55).  The
random batch is deterministic so its evidence can be reproduced exactly.
"""

from __future__ import annotations

import argparse
import json
import random
from math import gcd, prod
from pathlib import Path

import sympy as sp


N = 55
EXPECTED = 134694094094758395331307111329132
WINNER = "0000000100011111011011001101011011010100011110101000111"
COMPLEMENT = "".join("1" if bit == "0" else "0" for bit in WINNER)
FOURIER_PRIMES = (1000000321, 1000000871, 1000001311)


def is_prime_trial(value: int) -> bool:
    if value < 2:
        return False
    if value % 2 == 0:
        return value == 2
    divisor = 3
    while divisor * divisor <= value:
        if value % divisor == 0:
            return False
        divisor += 2
    return True


def circulant(word: str) -> list[list[int]]:
    if len(word) != N or set(word) - {"0", "1"}:
        raise ValueError("expected a 55-bit word")
    return [[int(word[(column - row) % N]) for column in range(N)]
            for row in range(N)]


def bareiss(matrix: list[list[int]]) -> int:
    """Fraction-free determinant with exact Python integers."""

    size = len(matrix)
    if size == 0:
        return 1
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for pivot_index in range(size - 1):
        pivot_row = next(
            (row for row in range(pivot_index, size)
             if a[row][pivot_index] != 0),
            None,
        )
        if pivot_row is None:
            return 0
        if pivot_row != pivot_index:
            a[pivot_index], a[pivot_row] = a[pivot_row], a[pivot_index]
            sign = -sign
        pivot = a[pivot_index][pivot_index]
        for row in range(pivot_index + 1, size):
            for column in range(pivot_index + 1, size):
                numerator = (
                    a[row][column] * pivot
                    - a[row][pivot_index] * a[pivot_index][column]
                )
                if numerator % previous:
                    raise ArithmeticError("Bareiss non-exact division")
                a[row][column] = numerator // previous
        for row in range(pivot_index + 1, size):
            a[row][pivot_index] = 0
        previous = pivot
    return sign * a[-1][-1]


def gaussian_mod(matrix: list[list[int]], modulus: int) -> int:
    """Determinant by row-reduction in a prime field."""

    a = [[entry % modulus for entry in row] for row in matrix]
    result = 1
    size = len(a)
    for column in range(size):
        pivot = next((row for row in range(column, size)
                      if a[row][column]), None)
        if pivot is None:
            return 0
        if pivot != column:
            a[column], a[pivot] = a[pivot], a[column]
            result = (-result) % modulus
        diagonal = a[column][column]
        result = result * diagonal % modulus
        inverse = pow(diagonal, modulus - 2, modulus)
        for row in range(column + 1, size):
            factor = a[row][column] * inverse % modulus
            if not factor:
                continue
            for tail in range(column + 1, size):
                a[row][tail] = (a[row][tail] - factor * a[column][tail]) % modulus
    return result


def primitive_root_55(modulus: int) -> int:
    if (modulus - 1) % N:
        raise ValueError("modulus is not 1 modulo 55")
    for generator in range(2, 1000):
        root = pow(generator, (modulus - 1) // N, modulus)
        if (pow(root, N, modulus) == 1
                and pow(root, N // 5, modulus) != 1
                and pow(root, N // 11, modulus) != 1):
            return root
    raise RuntimeError("no primitive 55th root found in search range")


def fourier_tables(modulus: int) -> list[list[int]]:
    root = primitive_root_55(modulus)
    return [[pow(root, frequency * position, modulus)
             for position in range(N)] for frequency in range(N)]


def fourier_mod(word: str, modulus: int, table: list[list[int]] | None = None) -> int:
    powers = table if table is not None else fourier_tables(modulus)
    values = [sum(int(bit) * powers[frequency][position]
                  for position, bit in enumerate(word)) % modulus
              for frequency in range(N)]
    return prod(values, start=1) % modulus


def resultant_product(word: str) -> tuple[tuple[int, int, int], int]:
    x = sp.Symbol("x")
    polynomial = sum(int(bit) * x**position for position, bit in enumerate(word))
    norms = tuple(
        int(sp.resultant(polynomial, sp.cyclotomic_poly(order, x), x))
        for order in (5, 11, 55)
    )
    return norms, sum(map(int, word)) * prod(norms)


def affine_orbit(word: str) -> set[str]:
    return {
        "".join(word[(unit * index + translation) % N]
                for index in range(N))
        for unit in range(1, N)
        if gcd(unit, N) == 1
        for translation in range(N)
    }


def check_one(word: str, tables: dict[int, list[list[int]]]) -> dict[str, object]:
    matrix = circulant(word)
    exact = bareiss(matrix)
    norms, resultant = resultant_product(word)
    modular = {
        str(modulus): gaussian_mod(matrix, modulus)
        for modulus in FOURIER_PRIMES
    }
    fourier = {
        str(modulus): fourier_mod(word, modulus, tables[modulus])
        for modulus in FOURIER_PRIMES
    }
    for modulus in FOURIER_PRIMES:
        key = str(modulus)
        expected_mod = exact % modulus
        if modular[key] != expected_mod or fourier[key] != expected_mod:
            raise AssertionError((word, modulus, exact, modular[key], fourier[key]))
    if resultant != exact:
        raise AssertionError((word, exact, norms, resultant))
    if exact < 0:
        raise AssertionError((word, "odd-order binary determinant must be nonnegative"))
    return {
        "word": word,
        "weight": word.count("1"),
        "bareiss": exact,
        "norms": list(norms),
        "resultant_product": resultant,
        "modular_gaussian": modular,
        "fourier_product": fourier,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--random-count", type=int, default=100)
    parser.add_argument("--output", type=Path,
                        default=Path(__file__).with_name("winner_audit.json"))
    args = parser.parse_args()
    if args.random_count < 100:
        raise SystemExit("--random-count must be at least 100")
    if any(not is_prime_trial(p) or p % 55 != 1 for p in FOURIER_PRIMES):
        raise AssertionError("independent primes failed primality/congruence check")
    tables = {p: fourier_tables(p) for p in FOURIER_PRIMES}
    if any(pow(tables[p][1][1], 55, p) != 1 for p in FOURIER_PRIMES):
        raise AssertionError("Fourier table root check failed")

    winner = check_one(WINNER, tables)
    complement = check_one(COMPLEMENT, tables)
    if winner["bareiss"] != EXPECTED:
        raise AssertionError(winner)
    if complement["bareiss"] * 28 != EXPECTED * 27:
        raise AssertionError("complement duality ratio failed")

    orbit = affine_orbit(WINNER)
    stabilizer = sum(
        "".join(WINNER[(unit * index + translation) % N]
                for index in range(N)) == WINNER
        for unit in range(1, N)
        if gcd(unit, N) == 1
        for translation in range(N)
    )
    if len(orbit) != 2200 or stabilizer != 1:
        raise AssertionError((len(orbit), stabilizer))
    orbit_residues = {
        p: {fourier_mod(word, p, tables[p]) for word in orbit}
        for p in FOURIER_PRIMES
    }
    if any(len(values) != 1 for values in orbit_residues.values()):
        raise AssertionError("affine orbit Fourier determinants disagree")

    rng = random.Random(550055)
    random_words = [format(rng.getrandbits(N), f"0{N}b")
                    for _ in range(args.random_count)]
    random_checks = [check_one(word, tables) for word in random_words]
    output = {
        "method": "clean-room Bareiss + modular Gaussian/Fourier + SymPy resultants",
        "primes": list(FOURIER_PRIMES),
        "winner": winner,
        "complement": complement,
        "complement_ratio": {
            "left": complement["bareiss"] * 28,
            "right": EXPECTED * 27,
        },
        "orbit_size": len(orbit),
        "stabilizer_size": stabilizer,
        "orbit_fourier_residue_cardinalities": {
            str(p): len(values) for p, values in orbit_residues.items()
        },
        "random_count": len(random_checks),
        "random_checks": random_checks,
    }
    args.output.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "winner": winner["bareiss"],
        "complement": complement["bareiss"],
        "random_count": len(random_checks),
        "orbit_size": len(orbit),
        "stabilizer_size": stabilizer,
        "status": "PASS",
    }, sort_keys=True))


if __name__ == "__main__":
    main()
