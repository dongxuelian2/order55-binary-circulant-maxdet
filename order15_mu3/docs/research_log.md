# Research log

## 2026-09-12

- Isolated the new investigation from the existing order-55 project.
- Directly checked the 2025 LAA paper and the 2023/2026 dissertation.
- Recovered the explicit `M15` exponent matrix from dissertation page 167.
- Implemented integer-only Eisenstein arithmetic, exact division,
  fraction-free determinant, exact Gram construction, dephasing, and the full
  row-inner-product catalogue.
- Added an exact benchmark verifier and regression tests. Verification
  reproduces `det=604661760+241864704 omega` and norm
  `277868041444786176`; the exact determinant of `HH*` agrees.
- The recovered matrix has 84 orthogonal row pairs, 14 nonzero Gram pairs of
  norm 3, and 7 of norm 9.
- Derived the rigorous upper bound `D_3(15) <= 3249*15^12` and the stronger
  conditional rule that a record improvement can contain no Gram entry of
  norm 84 or larger.
- Implemented and tested an inverse-update single-entry ascent engine. An
  initial 200-basin iterated search (seed 1503, perturbations of the dephased
  benchmark plus periodic random starts) found no improvement; this is a
  negative heuristic result, not evidence of maximality.
- A further 4,000 entry-only basins (seed 1515) and 500 hybrid basins using
  globally optimal row/column best responses (seed 1516) also found no
  improvement. The hybrid row optimizer was separately checked against brute
  force on randomized small instances.
- Scanned all `3^14=4,782,969` dephased circulant first rows by vectorized FFT
  ranking and exactly certified 128 finalists. The best finalist norm was
  `252876769144602624`, below the unrestricted record. Because candidate
  ranking used floating point, this run is not yet an exact proof of the
  circulant optimum, though it reproduces the scale of the thesis record.
- Proved by exact trace optimization that a strict improvement has total Gram
  energy `Q<=168`, and derived the row-color congruence distinguishing norms
  divisible by 9 from norms congruent to 3 modulo 9.
- Exactly excluded every counterexample with friendship support `F_7` and
  reduced the possible row-color partitions to `(15,0,0)`, `(14,1,0)`,
  `(13,2,0)`, and `(13,1,1)`.
- Located the exact coding-theory result `B_3(15,10)=12` (Todorov--Bogdanova,
  DOI `10.28919/jmcs/4964`) and converted it into internal-support
  independence constraints.  This excludes 14-class internal energies 0 and
  9, the energy-18 path, and the energy-27 one-edge and star profiles.
- Exhaustively enumerated abstract `(14,1,0)` Gram blocks through internal
  energy 27.  At energy 18 only two disjoint norm-9 edges survive; at energy
  27 only connected `P4` and triangle profiles retain norm-admissible values.
- For a 13-row class with minimum internal energy 9, exact two-hub Schur
  enumeration forces all 26 cross norms to be 3 and leaves 20 determinant
  values for `(13,2,0)` and 23 for `(13,1,1)`.
- Added PySAT encodings for partial orthogonal systems.  The encoding
  constructs nine pairwise orthogonal rows, proves that one such nine-row seed
  is maximal, and reduces the 13-row problem to five contingency-table orbits;
  one orbit is UNSAT and four remain timeout-unknown.
- Exhaustively classified the all-same-color support through energy 45.
  Energies 27 and 36 have no arithmetic survivor.  At energy 45 exactly eight
  determinant values remain, on `K3+2K2`, `C5`, or triangle-with-length-2-tail
  support.
- Derived the row/column Gram intertwining identity and an exact third-root
  kernel-label propagation algorithm.  This excludes every all-same-color
  case at energies 45, 54, 63, and 72; the last energy-54/72 `K4` and
  energy-63 five-vertex profiles force non-PSD residual Gram matrices.
- Completed the all-same-color energy-81 layer.  A two-sided nullity/rank
  inequality eliminates every isolated support.  The only no-isolate support
  is `K3+6K2`; spectral and phase-balance reductions leave one `3+6*2` block
  CNF, independently certified UNSAT by Glucose and Maple.  Hence this color
  partition has energy at least 90.
- Found and proved a twin-leaf PSD obstruction for `(14,1,0)`, excluding total
  energies 60 and 69.  Exactly enumerated all three total-energy-78 slices;
  their arithmetic maxima are below the trace bound at energy 87.
- Improved the rigorous global integer upper bound to
  `337136180404499901`, a factor `1.2132959899` above the published record.
- Corrected a SAT symmetry subtlety: a row fixed as a contingency-table orbit
  representative is no longer also forced to be lexicographically first.
  Pre-correction orbit-UNSAT claims were withdrawn; fixed-seed UNSAT and the
  new reduced energy-81 UNSAT certificate are unaffected.
- Completed exact all-same support certificates at Q=90, Q=99, and Q=108.
  Actual `mu3` kernel balance plus two-sided nullity/rank pairing eliminates
  every isolated support.  Complete no-isolate enumeration gives caps
  `288608006862471168` and `282438239494864896` at Q=90 and Q=99, and no
  Q=108 arithmetic value above the record.
- Added Schur/color certificates at Q=96, Q=105, and Q=114.  Tight sparse
  blocks use a clique/Motzkin--Straus spectral estimate; tight size-13 blocks
  retain the mandatory outside edge.  The rigorous upper is now
  `312967780240536567`, gap factor `1.1263180127273646`.
- Began Q=117: all mixed isolated and mixed no-isolate partitions are empty
  above the record.  The targeted no-isolate thirteen-unit enumeration used
  6,415 component types and 10,329 profiles and found no arithmetic survivor.
  The isolated thirteen-unit and size-13 Q=117 layers remain active.
