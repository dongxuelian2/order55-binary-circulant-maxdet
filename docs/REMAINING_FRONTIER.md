# Remaining mathematical frontier

The high-energy tail is now closed completely.  The deterministic certificates
exclude every genuine shell `Q=153,159,162,168`, while trace stability already
puts `Q>=171` below the record.  Consequently every strict counterexample must
satisfy

`Q <= 150`.

The post-tail aggregate has now been recomputed.  The exact record is

`277868041444786176 <= D_3(15)`

and the current rigorous aggregate upper is

`D_3(15) <= 287687522570016537`.

In the Hadamard-normalized scale
`R_15 = 100 |det H| / 15^(15/2)`, this is

`79.6590012604% <= R_15 <= 81.0543020059%`.

The old Q=153 trace envelope no longer controls the upper bound.  The new
numerical bottleneck is a Q=126 all-same-color weighted-K5 branch with twelve
support edges, seven isolates, and eigenvalue cap `287/10`.  Its exact rational
KKT upper is approximately `2.876875225936136e17`; after the mandatory `3^14`
divisibility and rational Eisenstein-norm sieve, four non-norm multiples are
skipped and the integer upper above results.

This creates an important distinction.  **Q=150 is the highest unresolved
energy shell**, so it remains the natural top-down structural frontier.
However, **Q=126 currently controls the global numerical upper**, so improving
that weighted-K5 certificate is the most direct way to narrow the displayed
upper/lower interval.

## Final tail closure

The row/column color congruence gives

- `(15,0,0)` and `(13,1,1)`: `Q ≡ 0 (mod 9)`;
- `(14,1,0)` and `(13,2,0)`: `Q ≡ 6 (mod 9)`;
- `Q ≡ 3 (mod 9)` is impossible.

Together with the trace-stability cutoff, the only genuine post-Q150 shells
were therefore `Q=153,159,162,168`.

`order15_mu3/scripts/verify_tail_structural_reduction.py` first closes Q=162
and Q=168 and reduces Q=153/Q=159 to nine explicit branches.
`order15_mu3/scripts/verify_tail_final_closure.py` closes those last branches.

### Sparse size-13 branches at Q=153 and Q=159

For a `13+2` Schur decomposition write

`S = [[15-x, g-z], [conj(g-z), 15-y]]`,

where `x=u* A^{-1}u`, `y=v* A^{-1}v`, and `z=u* A^{-1}v`.
Cauchy--Schwarz gives `|z|<=sqrt(xy)`.  Setting `s=x+y`, the determinant is
bounded by

`det(S) <= (15-s/2)^2 - max(0, |g|-s/2)^2`.

This expression decreases in the relevant range, while
`s >= cross_energy/lambda_max(A)`.  The exact sparse-block enumeration and
rational spectral brackets from the structural certificate therefore turn the
previous five residual size-13 allocations into strict below-record bounds.
This closes all `(13,1,1)` cases at Q=153 and the final `(13,2,0)` case at
Q=159.

### Q=159, color type (14,1)

The equal-color-sign orbit lemma excludes internal support isolates.  Split the
14-row internal support into connected components.  A component on `v`
vertices with `u` nine-energy units must satisfy `u>=v-1`.  Its determinant is
bounded by the exact stationary product in dimension `v`, while

`lambda_max <= 15 + sqrt(18 u (v-1)/v)`.

The singleton cross vector contributes norm at least three on every vertex, so
Schur correction can be bounded component by component.  Exhausting the finite
component allocations closes internal energies 99 and 108 directly.  At
internal energy 117, every disconnected allocation is also closed; the only
remaining abstract configuration is one connected 14-vertex component with 13
unit edges.  It is therefore a tree with every edge of norm 9.

There are exactly 3,159 unlabeled trees on 14 vertices.  Their Gram
determinants are evaluated exactly through the tree matching polynomial.  The
maximum is attained uniquely by the path `P14`, with determinant

`16802420983158456`,

and its Schur bound is still strictly below the record.  Thus Q=159 is closed.

### Q=153, all-same color

After the size-13 cases are closed, both row and column color partitions of a
Q=153 counterexample would have to be all-same.  Write

`G=15I+3A`, `K=15I+3B`, so `AH=HB`.

If the row support has at most four isolates, a K6 would consume 15 of the 17
available support-energy units and leave too little energy to cover the other
active vertices.  Hence the clique number is at most five.  Weighted
Motzkin--Straus plus the capped-eigenvalue KKT certificate places this case
below the record.

If there are at least eight isolates, the corresponding rows of the invertible
matrix H are independent vectors in `ker(B)`.  Thus K has at least eight
eigenvalues exactly 15.  Fixing those eigenvalues and optimizing the remaining
seven at the Q=153 variance is already below the record.

The intermediate isolate counts 5, 6, and 7 are the only K6 boundary cases.
The 17-unit budget forces the K6 itself to have fifteen norm-9 edges and leaves
only one or two unit edges outside it.  Fischer and two-dimensional Schur bounds
then close all three cases exactly.  Hence Q=153 is closed.

## New frontier

The strongest tail statement remains

`strict counterexample => Q <= 150`.

For structural closure, continue with Q=150 and descend only when necessary.
For the fastest improvement of the numerical global bound, sharpen the Q=126
twelve-edge weighted-K5 / seven-isolate branch first.  Both directions reuse
existing exact support, Schur, and capped-eigenvalue machinery.

All final claims remain exact integer/rational certificates; floating-point
spectra are not proof inputs.
