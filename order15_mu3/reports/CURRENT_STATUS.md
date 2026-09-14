# Current status (2026-09-14)

- **Current exact record / lower bound:** recovered and exactly verified published `M15`,
  `det=604661760+241864704 omega` and
  `|det|^2 = 2^22 3^20 19 = 277868041444786176`.
- **Post-tail rigorous aggregate upper:** `287687522570016537` for `|det H|^2`.
  In the Hadamard-normalized scale `100 |det H| / 15^(15/2)`, the rigorous
  interval is
  `79.6590012604% <= R_15 <= 81.0543020059%`.
- **Arithmetic post-processing:** the controlling strict rational bound is
  rounded down to the largest possible integer, restricted to multiples of
  `3^14`, then filtered by the rational Eisenstein-norm criterion.  Four
  divisible non-norm candidates are skipped before reaching the displayed
  integer upper.
- **Numerical bottleneck:** after removing the post-Q150 tail, the largest
  surviving analytic source is no longer Q=153 or Q=150.  It is a Q=126
  all-same-color weighted-K5 branch with twelve support edges, seven isolates,
  and the exact rational eigenvalue cap `287/10`.  Its KKT product is about
  `2.876875225936136e17`, slightly above the Q=150 boundary bound.
- **Energy frontier versus numerical bottleneck:** the highest energy at which
  a strict counterexample can still occur is Q=150, but the loosest currently
  certified determinant upper occurs at Q=126.  These are different notions.
- **Exact maximality proved:** no.
- **Strongest verified tail theorem:** every strict counterexample satisfies
  `Q<=150`.  The only genuine post-Q150 shells were `Q=153,159,162,168`; all
  four are now strictly below the published record, while trace stability
  already excludes `Q>=171`.
- **Q=162 and Q=168:** closed by
  `scripts/verify_tail_structural_reduction.py`, using color congruence, exact
  sparse size-13 abstract-Gram enumeration, rational Sylvester spectral
  brackets, Schur residual bounds, exact internal determinant maxima, and
  weighted Motzkin--Straus spectral caps.
- **Q=153 closure:** the former sparse `(13,1,1)` residual allocations are
  eliminated by a two-hub Schur--Cauchy inequality; the all-same branch is
  closed by the clique/isolate/intertwining split in
  `scripts/verify_tail_final_closure.py`.
- **Q=159 closure:** the final sparse `(13,2,0)` allocation is removed by the
  same Schur--Cauchy bound.  The `(14,1,0)` boundary reduces to all 3,159
  unlabeled order-14 trees; exact matching-polynomial enumeration has unique
  maximum at `P14`, still below the record.
- **Next research targets:** Q=150 remains the highest unresolved shell, while
  the Q=126 twelve-edge weighted-K5 branch is the best target for lowering the
  *global numerical upper bound*.  A final proof must eventually bring every
  remaining record-level branch down to the published record.
