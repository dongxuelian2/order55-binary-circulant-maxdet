# Tail closure certificate map — 2026-09-14

This note records the exact logical dependency of the completed post-Q150 tail
closure for the order-15 third-root maximal-determinant problem.

1. `verify_trace_stability.py` proves every strict counterexample has `Q<=168`
   and excludes `Q>=171`.
2. Color congruence leaves only `Q=153,159,162,168` above Q=150.
3. `verify_tail_structural_reduction.py` closes Q=162 and Q=168 and reduces
   Q=153/Q=159 to nine explicit branches.
4. `verify_tail_final_closure.py` applies the two-hub Schur--Cauchy bound to the
   five sparse size-13 residuals, closes Q=159 type `(14,1)` by componentwise
   Schur bounds plus exact enumeration of all 3,159 order-14 trees, and closes
   the Q=153 all-same branch by the clique/isolate/intertwining split.
5. Therefore every strict counterexample satisfies `Q<=150`.

The record remains

`|det(H)|^2 = 277868041444786176 = 2^22 * 3^20 * 19`.

This tail theorem does **not** by itself prove exact maximality.  The present
Q=150 certificate was constructed as a comparison against the former Q=153
trace envelope.  The next proof obligation is a direct record-level Q=150
closure, followed by recomputation of the aggregate arithmetic upper bound.
