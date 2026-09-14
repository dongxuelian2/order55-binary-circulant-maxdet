# Current status (2026-09-14)

- **Current exact record:** recovered and exactly verified published `M15`,
  `det=604661760+241864704 omega` and
  `|det|^2 = 2^22 3^20 19 = 277868041444786176`.
- **Rigorous upper bound:** `287710229992756239`, from the full unmocked
  pre-Q144 replay, exact Q=144/Q=150 branch certificates, the Q=153 trace
  envelope, `3^14` divisibility, and the Eisenstein norm sieve.
- **Exact maximality proved:** no.
- **Strongest verified tail theorem:** after the exact `Q<=150` branch
  certificates, the structural tail replay excludes the complete `Q=162` and
  `Q=168` shells.  Hence every strict counterexample in the high-energy tail
  has `Q in {153,159}`.  The proof combines the color congruence, exact sparse
  size-13 abstract-Gram enumeration, rational Sylvester spectral brackets,
  Schur residual bounds, exact internal determinant maxima, and weighted
  Motzkin--Straus spectral caps.  Replay:
  `scripts/verify_tail_structural_reduction.py`.
- **Q=153 residual frontier:** the all-same-color branch plus exactly four
  `(13,1,1)` allocations remain after the new certificate:
  `(e,c,r)=(27,78,48),(27,87,39),(36,78,39),(54,78,21)`.
- **Q=159 residual frontier:** exactly one `(13,2,0)` allocation remains,
  `(e,c,r)=(54,78,27)`, together with the `(14,1,0)` internal-energy slices
  `e=99,108,117`.  Thus the entire post-Q150 problem has been reduced to nine
  explicit branches: five at Q=153 and four at Q=159.
- **Q=162 closure:** the all-same-color support has clique number at most six;
  the resulting capped-eigenvalue KKT product is below the record.  Every
  `(13,1,1)` allocation is below the record after the exact two-hub Schur
  residual replay.
- **Q=168 closure:** every `(13,2,0)` allocation is below the record after the
  same size-13 replay, and every `(14,1,0)` allocation is below the record by
  the fourteen-dimensional internal spectral/Schur certificate.
- **Global numerical gap:** the upper/lower ratio remains
  `60153061831/58095304704 ≈ 1.03542036895`, with additive squared gap
  `9842188547970063`.  This number is unchanged because the surviving Q=153
  trace envelope still controls the aggregate upper bound.
- **Completed Q=117 partial searches:** every mixed isolated and mixed
  no-isolate partition is empty above the record; the no-isolate thirteen-unit
  layer has 10,329 exact profiles and no arithmetic survivor.  These remain
  useful supporting certificates for the lower-energy structure.
