# Current status (2026-09-14)

- **Current exact record:** recovered and exactly verified published `M15`,
  `det=604661760+241864704 omega` and
  `|det|^2 = 2^22 3^20 19 = 277868041444786176`.
- **Committed aggregate upper before this tail closure:**
  `287710229992756239`, from the pre-Q144 replay, exact Q=144/Q=150 boundary
  certificates, the former Q=153 trace envelope, `3^14` divisibility, and the
  Eisenstein norm sieve.  This aggregate number has not yet been recomputed
  against the new tail theorem because the next controlling frontier is Q=150.
- **Exact maximality proved:** no.  The tail is now finished, but the existing
  Q=150 certificate compares that shell with the old Q=153 envelope rather
  than directly with the record.
- **Strongest verified tail theorem:** every strict counterexample satisfies
  `Q<=150`.  The only genuine post-Q150 shells were `Q=153,159,162,168`; all
  four are now strictly below the published record, while trace stability
  already excludes `Q>=171`.
- **Q=162 and Q=168:** closed by
  `scripts/verify_tail_structural_reduction.py`, using color congruence, exact
  sparse size-13 abstract-Gram enumeration, rational Sylvester spectral
  brackets, Schur residual bounds, exact internal determinant maxima, and
  weighted Motzkin--Straus spectral caps.
- **Q=153 sparse size-13 closure:** the four former `(13,1,1)` residual
  allocations `(e,c,r)=(27,78,48),(27,87,39),(36,78,39),(54,78,21)` are
  eliminated by a two-hub Schur--Cauchy inequality retaining the mandatory
  outside Gram norm.
- **Q=153 all-same closure:** if the support has at most four isolates, the
  17-unit energy budget forbids K6 and the clique-five spectral cap is below
  the record.  At least eight isolates force at least eight 15-eigenvalues on
  the opposite Gram through `AH=HB`, again below the record.  The intermediate
  isolate counts 5, 6, and 7 reduce to exact K6 boundary supports and are
  closed by Fischer/Schur bounds.
- **Q=159 closure:** the final `(13,2,0)` sparse allocation `(54,78,27)` is
  removed by the same Schur--Cauchy bound.  For `(14,1,0)`, componentwise
  Schur bounds close internal energies 99 and 108 and every disconnected
  energy-117 allocation.  The sole connected boundary is a 14-vertex tree
  with thirteen norm-9 edges.  Exact matching-polynomial enumeration of all
  3,159 unlabeled trees gives the unique maximum at `P14`, with internal Gram
  determinant `16802420983158456`; its final Schur bound is below the record.
  Replay: `scripts/verify_tail_final_closure.py`.
- **New frontier:** Q=150.  The natural continuation is to reuse the new
  Schur--Cauchy/component methods on the Q=150 residual allocations and compare
  them directly with the published record.  If Q=150 closes, the aggregate
  certificate should then be recomputed and the next lower boundary attacked
  only if it controls the resulting upper.
- **Completed Q=117 partial searches:** every mixed isolated and mixed
  no-isolate partition is empty above the record; the no-isolate thirteen-unit
  layer has 10,329 exact profiles and no arithmetic survivor.  These remain
  useful supporting certificates for the lower-energy structure.
