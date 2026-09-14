# Proof roadmap

The archived argument is a chain of strict upper bounds.  Every arrow is
replayed by a script and a pytest module; the chain ends at a continuous trace
envelope, not at equality with the record.

1. **Algebraic normalization.**  Dephase a μ₃ matrix, construct `G=HH*`, and
   use exact Eisenstein arithmetic.  Mandatory `π^14` divisibility gives
   `3^14 | |det H|²`.
2. **Record and coarse restrictions.**  The recovered `M15` supplies the
   lower bound.  Pairwise Fischer bounds exclude Gram norms ≥84, and trace
   stability excludes `Q≥171`.
3. **Support and color reduction.**  The norm-mod-9 row-color congruence,
   friendship exclusion, and `B₃(15,10)=12` reduce strict counterexamples to
   four color partitions and sparse support components.
4. **Exact low/mid-energy certificates.**  SAT/SMT certificates, component
   catalogues, inverse-loss tables, spectral caps, and Schur complements close
   Q=60 through Q=150 in stages.  Unknown solver results remain explicitly
   unknown; only exact arithmetic/UNSAT runs enter the theorem chain.
5. **Aggregate.**  Take the maximum of the exact branch bounds, Q=144 and
   Q=150 boundaries, and the Q≥153 trace envelope.  Round down using
   `3^14` divisibility and reject non-Eisenstein norms.  This is a rigorous
   integer upper, not a candidate construction.
6. **Missing step for maximality.**  Close the genuine Q=153,159,162,168 tail
   by a structural realizability theorem, or find a sharper spectral/phase
   constraint.  The aggregate verifier also checks Q=156 and Q=165
   conservatively because its loop is generic; those shells are impossible by
   color congruence and do not affect the bound.  A
   successful result must be independently replayable and preserve the exact
   integer-norm sieve.
