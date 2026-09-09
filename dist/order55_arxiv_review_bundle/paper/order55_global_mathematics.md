# Exact reductions for odd order 55

This is a derivation log, not a completed global theorem.

## Fourier convention and complement

Use C(a)[i,j]=a[(j-i) mod n]. Its eigenvalues are f(zeta^j), up to
permuting j. For odd n, real coefficients pair all nonzero frequencies.
Thus det C(a)=k product_{j=1}^{(n-1)/2}|f(zeta^j)|^2 >=0.
For 0<k<n, complement multiplies each nonzero eigenvalue by -1 and the
zero eigenvalue by (n-k)/k. There are n-1 even such signs, so
k det C(1-a)=(n-k)det C(a). Complement is a bijection between weight
strata. At n=55 both constant words have zero determinant. Consequently
D01(55)=max_{1<=k<=27}(55-k)D_circ(55,k)/k.

Translation of a by t multiplies the Fourier eigenvalues by roots whose
product is zeta^{t n(n-1)/2}=1 for odd n. Multiplication of indices by a
unit permutes the frequencies. Therefore the full affine action preserves
the determinant. At 55 its order is 2200; orbit sizes must be divided by
actual stabilizers. Correlations ignore translations and also identify
opposite shifts, giving an effective unit group of order 20.

## First and second moments of squared magnitudes

Put q_j=|f(zeta^j)|^2 and c_s=sum_t a_t a_{t+s}.
Orthogonality gives 2 sum q_j=55k-k^2, 2 sum_{s=1}^{27} c_s=k(k-1),
and k^4+2 sum q_j^2=55(k^2+2 sum_{s=1}^{27}c_s^2).
Among h integers of fixed nonnegative sum T, write T=h b+r. Their minimum
square sum is (h-r)b^2+r(b+1)^2, by transferring one from any two entries
differing by at least two. These identities use no floating point.

For nonnegative x_i with sum S and sum of squares at least L, contract
toward their common mean until squares reach max(L,S^2/h). The product
increases: concavity of log proves this for positive entries, and continuity
handles zero entries. At a positive constrained extremum the Lagrange
equation 1/x=alpha+2 beta x has at most two roots. If the low root has
m occurrences, its deviation squared is V(h-m)/(hm); the high deviation
squared is Vm/[h(h-m)], V=L-S^2/h. Enumerating m=1,...,h-1 covers every
interior extremum. Boundary products are zero. Rational outward intervals
for both roots give a rigorous upper bound. The implementation is derived
from these equations; ALETHEIA-54 bounds.py uses the same extremal method.
Ryser's bound here also follows directly from the first-moment AM-GM:
det C(a)<=(k)[k(55-k)/54]^27. For the higher complementary weight replace
the first factor k by 55-k. Integer flooring is legitimate.

## Cyclotomic sectors and folds

The divisors 1,5,11,55 have totient dimensions 1,4,10,40, summing to 55.
The resultant product gives det C(a)=k N5 N11 N55. Each nontrivial
cyclotomic polynomial has even degree, so swapping the resultant operands
does not introduce a sign. All three norms are nonnegative.

For m=5 or 11, r_i=sum_{t=i mod m}a_t, and h_i=sum_{t=i mod m}c_t is
the circular autocorrelation of r. Orthogonality at order m gives
S1_m=(m sum r_i^2-k^2)/2 and S2_m=(m sum h_i^2-k^4)/2.
Subtract the 5 and 11 sectors from the full 55 moments for the primitive
55 sector. It has 20 pairs, while the other two have 2 and 5 pairs.
In particular, given both folds r,s its first moment is
T=(55k+k^2-5 sum r_i^2-11 sum s_j^2)/2, hence the complementary determinant
is at most (55-k) N5(r) N11(s) (T/20)^20.

Correlation sums over paired shifts divisible by 5 and by 11 are respectively
A=(sum r_i^2-k)/2 and B=(sum s_j^2-k)/2. They occupy 5 and 2 coordinates.
The remaining 20 have sum k(k-1)/2-A-B. Applying the integer square minimum
separately is a valid extra fourth-moment lower bound.

Every affine action on the mod-5 and mod-11 coordinates lifts jointly by
the Chinese remainder theorem. Enumerating all combinations of their
affine representatives therefore loses no determinant value. This does
NOT say that every pair of margin profiles has a binary lift.

## Why the (55,27,13) cyclic difference set cannot exist

A directly checkable algebraic proof avoids dependence on an existence table.
If c_s=13 for every nonzero s, put alpha=f(zeta_5). Then alpha*conjugate(alpha)=14.
Modulo 2, Phi5=x^4+x^3+x^2+x+1 is irreducible: a primitive fifth root in
characteristic 2 has degree ord_5(2)=4. Thus Z[zeta_5]/(2) is a field.
The displayed product is zero in that field, so alpha=0 modulo 2 (complex
conjugation is an automorphism). Hence alpha=2 beta with beta an algebraic
integer, forcing beta*conjugate(beta)=7/2. A rational algebraic integer is
an integer, contradiction. Thus E=sum(c_s-13)^2>0. Also E is even since
d_s^2=d_s mod 2 and sum d_s=0.

Published corroboration: S. Braic, "Primitive symmetric designs with at most
255 points", Glasnik Matematicki 45(65) (2010), pp.291-305, Theorem 4.4 and
p.297 explicitly discuss this parameter set. Apply p=2,w=5,j=2,n=k-lambda=14.
The theorem cites E. S. Lander, Symmetric Designs: An Algebraic Approach
(1983), pp.131,134. Primary PDF:
https://web.math.pmf.unizg.hr/glasnik/45.2/45(2)-01.pdf

## Correlation profiles and exact arithmetic

Even for integer profiles not known to be positive semidefinite, all q_j
are real. With known sum of their squares S2, product |q_j| is at most
(S2/27)^(27/2). This supplies a signed CRT capacity bound for the profile
enumerator independent of PSD or realizability. The configuration generator
checks the squared integer inequality against the product of the two primes.
Profiles passing the screen still require an exhaustive binary-lift proof.

## Final capped-moment refinement

The completed certificate additionally constrains every q_j by a uniform
rational cap from rearranging each fixed correlation multiset against the
three cosine coefficient multisets. To maximize a product on [0,B]^h,
enumerate the number of coordinates on the upper boundary B; the remaining
interior stationary coordinates again have at most two values. All zero
boundary products vanish. This cap excludes whole multisets which the
uncapped fourth moment leaves unresolved. See the final manuscript for the
proof and actual coverage counts. The final pi enclosure is obtained by
Machin's identity, not assumed decimal digits.
