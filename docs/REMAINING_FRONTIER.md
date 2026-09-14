# Remaining mathematical frontier

The low and middle energy slices are covered by exact support, Schur, SAT/SMT,
or rational spectral certificates through Q=150.  The high-energy continuation
now also closes the complete Q=162 and Q=168 shells.  Exact maximality remains
open, but the tail has been reduced from four energy shells to nine explicit
structural branches living only at Q=153 and Q=159.

## Tail-shell reduction

The proved row/column color congruence gives

- `(15,0,0)` and `(13,1,1)`: `Q ≡ 0 (mod 9)`;
- `(14,1,0)` and `(13,2,0)`: `Q ≡ 6 (mod 9)`;
- `Q ≡ 3 (mod 9)` is impossible.

Trace stability already gives `Q<=168` for a strict counterexample and puts
`Q>=171` below the published record.  Thus the original genuine tail shells
were `Q=153,159,162,168`.

The deterministic replay `order15_mu3/scripts/verify_tail_structural_reduction.py`
adds the following exact reductions.

### Q=162 is closed

For the all-same-color case, Q=162 contains at most eighteen unit support
weights.  A K7 would already require twenty-one unit edges, so the support
clique number is at most six.  Weighted Motzkin--Straus bounds the spectral
radius; the resulting capped-eigenvalue KKT determinant is strictly below the
record.

For `(13,1,1)`, write `Q=e+c+r`, where `e` is the internal energy of the
13-row class, `c` is the energy of its 26 cross entries, and `r` is the norm
between the two outside rows.  Sparse internal Gram blocks are enumerated
exactly.  Their smallest and largest eigenvalues are bracketed with rational
Sylvester tests, and the two-dimensional Schur complement retains the
mandatory outside residual.  Every admissible allocation lies below the
record.  Hence Q=162 cannot contain a strict counterexample.

### Q=168 is closed

The same exact size-13 replay closes every `(13,2,0)` allocation.  For
`(14,1,0)`, a fourteen-dimensional internal stationary determinant bound and
an exact rational spectral cap make every internal/cross allocation strictly
smaller than the record.  Hence Q=168 cannot contain a strict counterexample.

## Surviving Q=153 branches

Five branches remain:

1. the all-same-color `(15,0,0)` branch;
2. `(13,1,1)` with `(e,c,r)=(27,78,48)`;
3. `(13,1,1)` with `(27,87,39)`;
4. `(13,1,1)` with `(36,78,39)`;
5. `(13,1,1)` with `(54,78,21)`.

The all-same branch is already very close to the record under the capped
spectral product, but the current bound is not yet strict.

## Surviving Q=159 branches

Four branches remain:

1. `(13,2,0)` with `(e,c,r)=(54,78,27)`;
2. `(14,1,0)` with internal energy `e=99`;
3. `(14,1,0)` with internal energy `e=108`;
4. `(14,1,0)` with internal energy `e=117`.

All other Q=159 allocations are strictly below the record under the structural
certificate.

## What a successful continuation must do

The remaining work is now finite and sharply localized.  For the four sparse
size-13 allocations at Q=153 and the single one at Q=159, the natural next
step is to replace the spectral Schur relaxation by exact cross-vector phase
enumeration and then enforce the opposite Gram through `AH=HB`.  For the
three Q=159 size-14 slices, the internal no-isolate support can be enumerated
at energies 99, 108, and 117 and paired with the singleton cross vector.  The
Q=153 all-same branch should be attacked from both sides simultaneously using
support components, kernel/nullity data, and the intertwining rank inequality.

Any final closure should continue to emit exact integer/rational certificates;
floating-point eigenvalues may be used only for discovery.  The aggregate
integer upper remains `287710229992756239` until the Q=153 envelope itself is
strictly improved or eliminated.
