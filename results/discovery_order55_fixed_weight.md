# Discovery log: order-55 fixed-weight binary circulants

Date: 2026-09-09.  Host: AMD Ryzen 7 7840H, 8 physical / 16 logical cores.
Compiler: clang 22.1.7, `-O3 -std=c++17 -pthread`.  Search used 15 workers.

## Unrestricted record attacks

The type-1/type-3 bordered two-circulant objective was evaluated by Fourier
block diagonalization in O(m) per sign flip.  The executed algorithms were
record-seeded multi-start search, one-flip steepest quench, simulated annealing,
and iterated local search with variable-size kicks.

| order | type | seed | valid evaluations | wall seconds | best outcome |
|---:|---:|---:|---:|---:|---|
| 103 | 1 | 103002 | 430088932 | 90.0045 | exactly the public record |
| 111 | 3 | 111003 | 258413992 | 60.0032 | exactly the public record |
| 115 | 1 | 115004 | 247228002 | 60.0066 | exactly the public record |
| 119 | 3 | 119005 | 233913836 | 60.0058 | exactly the public record |
| total | | | 1169644762 | 270.0201 | no unrestricted improvement |

An earlier 45-second pilot at order 103 exposed an incorrect hand expansion
of the zero-frequency 3-by-3 determinant.  Its 254008430 scores are excluded
from every valid count above.  The correction was guarded by direct NumPy
matrix tests at public orders 103 and 111 before production was restarted.

All four valid runs returned the supplied record without a single accepted
final improvement.  This repeated ceiling across three block sizes and both
supported layouts triggered the predeclared switch to the exact structured
frontier.

## Exact structured exhaustion

For fixed weight k, every cyclic orbit of binary supports has a rotation that
contains position zero.  The enumerator therefore visits all `C(54,k-1)` such
supports.  It computes each determinant modulo the two primes
2305843009213696591 and 2305843009213697141 via the 55 Fourier eigenvalues in
the corresponding finite fields.  Signed CRT is unique because Hadamard gives
`|det|^2 <= k^55`, while half the CRT modulus is larger than `7^(55/2)`.

| weight | representatives | exact maximum | maximizing first row |
|---:|---:|---:|---|
| 1 | 1 | 1 | `1000000000000000000000000000000000000000000000000000000` |
| 2 | 54 | 2048 | `1000000000010000000000000000000000000000000000000000000` |
| 3 | 1431 | 1564031349 | `1000010000000001000000000000000000000000000000000000000` |
| 4 | 24804 | 6424318816256 | `1000010000010000000000000001000000000000000000000000000` |
| 5 | 316251 | 12556202952818275 | `1000000011000000000000000000000100000000000000001000000` |
| 6 | 3162510 | 4342461217049330166 | `1000000000001100100000100000000000000000000000001000000` |
| 7 | 25827165 | 570999857161750276477 | `1000000000000011000010000001000000001000000010000000000` |

Total exact representatives: 29332216.  The largest run (weight 7) took 8.8194
seconds at 2.93 million exact candidates per second.  `verify_structured.py`
independently reconstructed all seven matrices with SymPy, recomputed the exact
determinants, checked both residues and all combination counts, and returned
PASS for every row.
