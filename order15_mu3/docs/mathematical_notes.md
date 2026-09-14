# Mathematical notes

## Mandatory divisibility

Let `pi = 1 - omega`, so `N(pi)=3`. Every third root is congruent to `1`
modulo `pi`. Subtracting the first row of an order-`n` matrix from each of the
other `n-1` rows makes every entry of those rows divisible by `pi`. Hence
`pi^(n-1)` divides `det(H)` and `3^(n-1)` divides `N(det(H))=|det(H)|^2`.
At order 15 this proves the natural normalization by `3^14`.

The exponent is sharp as a universal divisibility statement. After dephasing,
take a core with `omega` on its diagonal and `1` elsewhere. Row subtraction
leaves the minor `(omega-1) I_(n-1)`, so the determinant has exactly
`pi^(n-1)` as its `pi`-part.

## Global obstruction to BH(15,3)

Equality in Hadamard's inequality would give `N(det H)=15^15=3^15 5^15`.
The rational prime 5 is inert in `Z[omega]` because `5 = 2 mod 3`; every
Eisenstein norm has even valuation at an inert prime. The odd valuation 15 is
impossible. This is global: individual row pairs can still be orthogonal,
because five copies of each root sum to zero.

## Order-15 row inner products

If the quotient counts are `(a,b,c)`, with `a+b+c=15`, the inner product is
`(a-c) + (b-c) omega` and its squared modulus is
`a^2+b^2+c^2-ab-ac-bc`. Exact enumeration gives 136 values. The least norms
are 0 (only `(5,5,5)`), then 3, 9, 12, and 21. The executable catalogue is
produced by `maxdet.mu3.inner_product_values(15)`.

## A first strict numerical upper bound

Let the support graph of `G-15I` join nonorthogonal row pairs and set
`r=sum_C(|C|-1)` over its nontrivial connected components. Order a spanning
forest parent before child and apply Gram--Schmidt. Each root contributes at
most 15; each nonroot has residual squared length at most
`15-3/15=74/5`, because every nonzero Gram entry has norm at least 3. Thus

`det(G) <= 15^(15-r) (74/5)^r`.

The cases `r=0` and `r=1` are impossible by Eisenstein norm parity at the
inert prime 5. Two disjoint edges (`r=2`) are impossible for the same reason:
their determinant is `15^11 (225-q1)(225-q2)`, and every allowed nonsingular
`225-q` has even 5-adic valuation.

The remaining `r=2` graph is one connected three-vertex component. If it is a
path, its determinant factor is `B=3375-15(q1+q2)`. For sums below 30 the only
possibilities are 6, 12, 15, 18, 21, 24; the corresponding factors
3285, 3195, 3150, 3105, 3060, 3015 each have an inert prime to odd exponent,
so none is a rational Eisenstein norm. Hence the path has `B<=2925`.
For a triangle put `s=q1+q2+q3`. AM--GM gives
`sqrt(q1 q2 q3)<=(s/3)^(3/2)`, so its determinant factor is at most
`3375-15s+2(s/3)^(3/2)`. This decreases on `9<=s<=675`, and is therefore
less than 3251. The factor is an integer norm; 3250 is not a norm
(`3250=2*5^3*13`), hence it is at most 3249. Finally for `r>=3`, the forest
bound is at most `222^3 15^9 < 3249 15^12`. Therefore

`D_3(15) <= 3249 * 15^12 = 57^2 * 15^12`.

This is rigorous but far above the benchmark and does not prove maximality.

## Immediate pairwise pruning above the benchmark

Fischer's inequality for a selected row pair gives
`det G <= (225-q)15^13`. With the verified benchmark, beating it requires
`q < 82.225`. Because allowed norms are discrete, every off-diagonal Gram
entry of any improvement must have norm at most 81; all norms 84 and above
are excluded.

## Trace-stability energy restriction

Put `Q=sum_(i<j)|G_ij|^2` and `S=tr((G-15I)^2)=2Q`. At a positive interior
maximum of the eigenvalue product with fixed trace and `S`, Lagrange
multipliers show that the eigenvalues take at most two values. If the high
value has multiplicity `m`, put `t^2=S/(15m(15-m))`. The stationary product
is `(15+(15-m)t)^m (15-mt)^(15-m)`. For each `m` it decreases with `t>0`;
boundary points have zero determinant.

Exact rational lower bounds on all 14 square roots at `S=342` certify that
every stationary product is strictly below the record. Since every allowed
Gram norm is a multiple of three, a hypothetical strict improvement must obey

`sum_(i<j)|G_ij|^2 <= 168`.

The rational replay is `scripts/verify_trace_stability.py`.

## Friendship support and a sign-orbit congruence

The verified record's support graph is the friendship graph `F_7`: one hub
joined to 14 leaves, plus a perfect matching on the leaves. Its 14 hub edges
have norm 3 and its seven matching edges have norm 9. After normalizing the
hub row to zero, any norm-3 leaf has exponent counts in one of two cyclic
sign orbits, distinguished by the sum of its 15 exponents modulo 3.

If two leaves are orthogonal, their difference counts are `(5,5,5)`, so their
exponent sums agree modulo 3. The cross-pair leaf graph of `F_7` is connected;
hence an all-norm-3 `F_7` profile forces all 14 leaves into the same sign
orbit. Global conjugation and individual leaf phase rotations then fix every
hub count triple to the same canonical representative. This is now imposed in
both exact realizability encodings.

In fact, for quotient counts `(a,b,c)` and row-color difference `delta`,

`delta = b+2c (mod 3)`, and `|G_ij|^2/3 = delta^2 (mod 3)`.

Hence equal-colored pairs have Gram norm divisible by 9, while
different-colored pairs have norm congruent to 3 modulo 9.

## Exact exclusion of friendship-support counterexamples

The congruence above completes the abstract `F_7` analysis. Cross-pair
orthogonality makes all leaves the same color. If the hub had that color too,
all 21 support edges would have norm at least 9, giving total energy at least
189, contrary to the counterexample bound 168. Thus hub-to-leaf edges have
norm `3 mod 9`, while matching edges have norm divisible by 9.

For each matched leaf pair, write its `2x2` block determinant as `d` and its
contribution to the hub Schur complement as `numerator/d`. Exhausting all
allowed Gram phases and norms at most 81 gives 872 local types, reduced by
exact dominance to

`(energy,d,numerator) = (15,216,72), (33,198,63), (42,189,54),`
`(69,162,45), (87,144,36)`.

Exact enumeration of all 330 multisets of seven frontier types shows that the
largest determinant is precisely `2^22 3^20 19`, attained only at the profile
using seven copies of `(15,216,72)`. Therefore a strict counterexample cannot
have nonorthogonality support graph `F_7`. The replay is
`scripts/verify_friendship_bound.py`.

## Global row-color partition reduction

If a color class of size `a` is internally orthogonal, split its Gram block
as `15 I_a` and take the Schur complement. Cross-color entries have norm at
least 3, so the Schur complement on the other `15-a` rows has trace at most
`(15-a)(15-a/5)`. AM--GM gives

`det G <= 15^a (15-a/5)^(15-a)`.

This is below the record for every `3<=a<=12`. Consequently every color
class in that size range must contain a same-color nonorthogonal pair, costing
at least 9 units of the 168 energy budget. Enumerating partitions and applying
refined trace/Schur bounds to the remaining size-10, size-11, and size-12
classes leaves only four possible color partitions for a strict counterexample:

`(15,0,0), (14,1,0), (13,2,0), (13,1,1)`.

All square-root comparisons and Schur inequalities are replayed with rational
arithmetic by `scripts/verify_color_partition_bounds.py`.

## Orthogonal-code bound and low-energy support classification

Pairwise orthogonal ternary phase rows differ in five coordinates by each of
the three phase quotients, and hence form a ternary equidistant Hamming code
of length 15 and distance 10.  Todorov and Bogdanova, *Ternary equidistant
codes of length 11<=n<=15* (2020), DOI `10.28919/jmcs/4964`, give the exact
value `B_3(15,10)=12`.  Therefore no 13 rows can be pairwise orthogonal.

Equivalently, the internal nonorthogonality support of a color class has
independence number at most 12.  Thus classes of sizes 13, 14, and 15 have
internal energy at least 9, 18, and 27 respectively.  For a 14-row class this
also classifies the smallest supports: energy 18 can only be two disjoint
norm-9 edges, while at energy 27 a single norm-27 edge and the three-edge star
are impossible.  The elementary graph replay is
`scripts/verify_support_independence.py`.

For color type `(14,1,0)`, exact block enumeration gives no arithmetic
candidate at internal energy 0, excludes energy 9 by the code bound, and at
energy 18 leaves only the two-disjoint-edge determinant list.  At energy 27,
the disconnected profiles have no Eisenstein-norm candidate; after the
independence reduction only the `P4` and triangle profiles retain candidates.

For either color type with a 13-row class, minimum internal energy is 9.  In
that case the class has one norm-9 edge and eleven isolated vertices.  Exact
Schur-loss enumeration shows that any cross norm greater than 3 among the 26
large-to-small pairs forces the determinant below the record.  With all 26
cross norms fixed at 3, oriented color-phase enumeration reduces the abstract
above-record lists to 20 rational Eisenstein norms for `(13,2,0)` and 23 for
`(13,1,1)`.  Replay: `scripts/verify_color_13_energy9.py`.

For the all-same-color partition, every energy is a multiple of 9.  Exact
unlabeled support-component enumeration first excludes energies 27 and 36.
At energy 45 arithmetic alone leaves eight values on `K3 + 2 K2`, `C5`, or a
triangle with a length-2 tail.

The decisive additional identity is obtained by writing

`G=HH*=15I+3A`, `K=H*H=15I+3B`.

Whenever both row and column partitions are all-same-color,
`A H = H B`.  An isolated row-support vertex therefore gives a full-support
third-root vector in the kernel of `B` on every nontrivial component.  Exact
unit-edge zero-sum enumeration plus this intertwining obstruction excludes
all energy-45 phases.  It reduces energy 54 to a singular weighted `K4` plus
11 isolates and energy 63 to one singular five-vertex graph plus 10 isolates.
In both cases the active kernel is one-dimensional.  The isolated rows then
force respectively `15I_11-4J_11` and `15I_10-5J_10`, with negative
all-ones eigenvalues `-29` and `-35`.  At energy 72 the unique no-isolate
support `P3+6K2` fails the Eisenstein norm test; the remaining last `K4` case
again has eigenvalue `-29`.  Thus energies 27 through 72 are all excluded.

At energy 81, isolated supports are eliminated by the two-sided rank bound

`m_row * v_col <= 15 d_col`,

where `m` is the number of isolates, `v` the number of active vertices, and
`d` the exact active-adjacency nullity.  The sole no-isolate graph is
`K3+6K2`.  Three of its four triangle phases fail spectral or exact
three-column balance constraints.  In the remaining negative-triangle phase,
intertwining forces a `3+6*2` block form: the triangle corner is constant,
outer triangle restrictions are permutations of the three roots, and every
outer `2x2` block is `[[x,y],[y,x]]`.  The complete reduced CNF is UNSAT.
Replay: `scripts/verify_color_15_low_energy.py`,
`scripts/verify_color_15_energy81.py`, and
`scripts/sat_q81_reduced.py`.  Consequently an all-same-color counterexample
has energy at least 90.

## Twin-leaf obstruction and refined global upper bound

For `(14,1,0)` at total energies 60 and 69, the minimum internal supports
have at least nine equal hub-only leaves after phase normalization.  If two
such leaf rows are subtracted, `(G-15I)H=H(K-15I)` places their difference in
the opposite support kernel.  The internal active blocks (two disjoint edges,
`P4`, or a triangle) are nonsingular, so the leaf rows agree on at least four
active columns.  Their complementary Gram would be `15I_m-vJ_m` with
`mv>15`, contradicting positive semidefiniteness.  Replay:
`scripts/verify_color_14_1_twin_bound.py`.

At total energy 78 the three possible splits `(internal,cross)` are
`(18,60)`, `(27,51)`, and `(36,42)`.  Exact arithmetic enumeration bounds
their largest norm-admissible determinants by, respectively,
`298188223593750000`, `307903629375000000`, and
`308697676962890625`.  All other surviving color cases have energy at least
87.  Continuing exact support and Schur certificates through Q=114, then
applying trace optimization from Q=117 through 168, followed by
`3^14` divisibility and the rational Eisenstein norm test, gives

`D_3(15) <= 312967780240536567`.

This is a squared-determinant gap factor `1.1263180127273646`.  At all-same
energies 90, 99, and 108, exact component enumeration and the actual `mu3`
kernel/rank obstruction eliminate isolated supports; complete no-isolate
enumerations place Q=90 and Q=99 below the next envelope and exclude every
Q=108 value above the record.  Sparse clique/spectral bounds and exact
outside-edge Schur estimates handle Q=96, 105, and 114.  Replay:
`scripts/verify_color_15_energy90.py`,
`scripts/verify_color_15_energy99.py`,
`scripts/verify_color_15_energy108.py`,
`scripts/verify_q96_boundary.py`, `scripts/verify_q105_boundary.py`,
`scripts/verify_q114_boundary.py`, and `scripts/verify_refined_upper.py`.
