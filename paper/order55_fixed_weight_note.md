# Exact determinants of sparse binary circulants of order 55

## Abstract

We determine the maximum absolute determinant of a 55-by-55 binary circulant
matrix at each fixed row weight from one through seven.  The computation
exhausts 29,332,216 rotation-covering representatives and evaluates every
determinant exactly with two finite-field Fourier transforms and signed Chinese
remainder reconstruction.  An independent SymPy verifier confirms each
maximizer.  The weight-seven maximum is
570999857161750276477.

## 1. Introduction and benchmark

For a binary word `a` of length n, write `circ(a)` for the circulant whose
rows are cyclic shifts of `a`, and define

`D_circ(n,k) = max |det circ(a)|`, over words of Hamming weight k.

Brent and Yedidia computed the unrestricted binary-circulant maximum through
order 53.  The ALETHEIA 54 project subsequently supplied an exact computational
proof at order 54.  A literature and data search completed on 2026-09-09 found
no published fixed-weight table at order 55 and found none of the integers in
the theorem below.  The scope here is deliberately structured: this note makes
no claim about the unrestricted 55-by-55 binary-circulant maximum.

## 2. Result

**Theorem.** For order 55 and weights one through seven,

| k | D_circ(55,k) |
|---:|---:|
| 1 | 1 |
| 2 | 2048 |
| 3 | 1564031349 |
| 4 | 6424318816256 |
| 5 | 12556202952818275 |
| 6 | 4342461217049330166 |
| 7 | 570999857161750276477 |

Maximizing first rows are stored in the machine-readable certificate.

## 3. Exhaustive method

Rotation of the first row changes a circulant determinant by at most sign.
Every nonempty support has a rotation containing coordinate zero.  Thus, for
weight k, it is sufficient (though not duplicate-free) to enumerate every
support containing zero: exactly `binom(54,k-1)` candidates.

Let `f(x)` be the binary polynomial associated with a support.  If a prime p
is 1 modulo 55 and omega has exact order 55 in the finite field, then

`det circ(a) = product_{j=0}^{54} f(omega^j) (mod p)`.

The implementation uses the independently primality-checked moduli
2305843009213696591 and 2305843009213697141.  It precomputes every required
root power and evaluates both residues for every enumerated support.  A signed
CRT value is unique: Hadamard's inequality gives
`|det circ(a)|^2 <= k^55`, the largest right-hand side occurs at k=7, and half
the CRT modulus has more than 121 bits of capacity versus a 77.21-bit bound.
Consequently the comparison of every reconstructed integer proves maximality
within each stated fixed-weight family.

## 4. Exact certificate and independent verification

The certificate records both primes, the two winning residues, the exact
representative count, each maximizing word, and each determinant.  The
independent verifier uses only the Python standard library and SymPy.  It
checks primality, congruence conditions, combination counts, weights, the
Hadamard/CRT inequality, reconstructs every 55-by-55 integer matrix, and
computes each determinant by exact domain elimination.  Its final output is
`order=55 total_representatives=29332216 PASS`.

To reproduce:

```powershell
clang++ -O3 -std=c++17 -Wall -Wextra -Wpedantic -pthread native/circulant_fixed_weight_exact.cpp -o native/circulant_fixed_weight_exact.exe
1..7 | ForEach-Object { native/circulant_fixed_weight_exact.exe --weight $_ --threads 15 }
.venv/Scripts/python.exe scripts/verify_structured.py
```

## 5. Reproducibility and limitations

The exhaustive generator is `native/circulant_fixed_weight_exact.cpp`; the
certificate is `certificates/order55_fixed_weight_1_7.json`.  The computation
ran on an AMD Ryzen 7 7840H using clang 22.1.7.  The seven strata contain
29,332,216 evaluated representatives in total; the largest stratum completed
in 8.8194 seconds with 15 workers.

This is an exact structured theorem, not a new unrestricted Hadamard maximal
determinant record and not the unrestricted binary-circulant value at order 55.
The novelty statement is conditioned on the documented public-source search;
absence from web-searchable sources is not a formal proof of nonpublication.

## References

1. R. P. Brent and A. B. Yedidia, “Computation of Maximal Determinants of Binary Circulant Matrices,” Journal of Integer Sequences 21 (2018), Article 18.5.6, https://arxiv.org/abs/1801.00399.
2. K. Terauchi, “ALETHEIA 54: exact discovery in binary circulants,” archived release, https://doi.org/10.5281/zenodo.21482404.
3. P. Browne, R. Egan, F. Hegarty, and P. Ó Catháin, “A Survey of the Hadamard Maximal Determinant Problem,” Electronic Journal of Combinatorics 28 (2021), P4.41, https://www.combinatorics.org/ojs/index.php/eljc/article/view/v28i4p41.
4. G. Butbaia et al., “New Records for the Hadamard Maximal Determinant Problem in Dimensions 51, 107, 111, 115, and 119,” arXiv:2608.22518v2, https://arxiv.org/abs/2608.22518.
