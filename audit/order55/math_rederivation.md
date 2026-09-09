# Independent mathematical rederivation

This note rederives the inequalities and quotient maps used by the order-55
certificate.  It is deliberately phrased independently of the production
Python modules.

## Circulant Fourier factorization

Let `a=(a_0,...,a_54)` be a binary word, let `A` be its 55 by 55 circulant,
and put

\[
k=\sum_t a_t,\qquad C_s=\sum_t a_ta_{t+s}\quad (s\bmod 55).
\]

Thus `C_0=k`, `C_{55-s}=C_s`, and
\[
\sum_{s=1}^{27}C_s=\frac{k(k-1)}2.
\]

For a primitive 55th root `omega`, the circulant eigenvalues are
\[
\lambda_j=\sum_{t=0}^{54}a_t\omega^{jt}.
\]
The zero frequency is `lambda_0=k`, while conjugacy gives
`lambda_{55-j}=conjugate(lambda_j)`.  Therefore
\[
\det A=k\prod_{j=1}^{27}|\lambda_j|^2\geq 0.
\]

The squared nonzero eigenvalues are determined by the half-correlation vector:
\[
|\lambda_j|^2
=k+2\sum_{s=1}^{27}C_s\cos(2\pi js/55).
\]
This makes the determinant a function of the correlation profile, before any
binary realizability assumption is made.

## Complement symmetry

For the complement `b=1-a`, every nonzero Fourier coefficient changes sign,
whereas the zero coefficient changes from `k` to `55-k`.  There are 54
nonzero frequencies, an even number, so
\[
k\det C(b)=(55-k)\det C(a).
\]
Consequently it is enough to control weights `0<=k<=27`; the endpoints have
zero determinant and the other weights are recovered by complementing.

## First and fourth moments

The mean of the 27 half-correlations is
\[
b_k=\frac{k(k-1)}{54}.
\]
Parseval at second order gives
\[
\sum_{j=1}^{27}|\lambda_j|^2=\frac{k(55-k)}2,
\]
so the first-moment AM-GM bound is
\[
\det A\leq k\left(\frac{k(55-k)}{54}\right)^{27}.
\]
The discrete Fourier Parseval identity applied to `C_s` gives
\[
\sum_{j=1}^{27}|\lambda_j|^4
=\frac{55\left(k^2+2\sum_{s=1}^{27}C_s^2\right)-k^4}{2}.
\]
If this sum is `S`, AM-GM gives the uncapped bound
\[
\det A\leq k\left(\frac{S}{27}\right)^{27}.
\]
The production certificate combines this with the first-moment bound and
uses the smaller outward-rounded integer floor.  Recomputing the rational
formula independently leaves only weights 25, 26, and 27 active.  Their
maximum allowed half-square sums are respectively 3338, 3919, and 4570.

## Spectral caps and exact rounding

For a sorted profile `c=(c_1,...,c_27)`, each eigenvalue square is bounded
using an upper/lower rational interval for
`2 cos(2*pi*r/55)`.  The audit recomputed pi with Machin's identity and a
longer alternating Taylor enclosure, then checked every stored dyadic
interval.  Products of rational interval endpoints were rounded outwards;
no floating-point comparison is used to discard a possible maximizer.

The same calculation is applied after each unit action on the shifts.  The
lexicographically least unit image is the canonical profile representative.
Units modulo 55 act through the 20 effective classes on paired shifts.

## Mod-5 and mod-11 folds

For `m` equal to 5 or 11, define the margin vector
\[
r_u=\sum_{t\equiv u\pmod m}a_t.
\]
Its entries sum to `k`, and its cyclic autocorrelation is the corresponding
fold of `C_s`.  The cyclotomic factorization
\[
\Phi_5\Phi_{11}\Phi_{55}=\frac{x^{55}-1}{x-1}
\]
splits the nontrivial Fourier product into the 5-, 11-, and 55-parts.  Exact
cyclotomic resultants of the margin polynomials provide the fold filters.
The audit independently enumerated the capped margin vectors, applied the
affine small-modulus canonicalization, and recomputed the signed two-prime
products.  The stored folded survivor sets and folded-correlation signatures
agree exactly.

## CRT exactness

The two order-55 fields use primes
`2305843009213696591` and `2305843009213697141`, both congruent to 1 modulo
55.  Independent modular arithmetic found primitive 55th roots in each
field.  The product of the primes is larger than the Hadamard bound for the
formal determinants being reconstructed, so the centered CRT representative
is the exact signed integer, not merely a residue.  The winner audit also
checks the same logic with three independent smaller primes and exact
integer/resultant calculations.

## Affine equivalence

Translations and multiplication by a unit modulo 55 permute rows and columns
of a circulant incidence matrix, preserving the determinant.  The affine group
therefore has size
\[
55\,\varphi(55)=55\cdot40=2200.
\]
The canonical winner has stabilizer 1 under this action, so its orbit contains
exactly 2200 words.  The independent orbit enumeration confirms closure and
the claimed size.

The remaining exhaustive step is not an assumption about profiles: each exact
55-profile target is joined with every compatible binary lift.  The independent
split-6 lift replay checks that final fibre search directly.
