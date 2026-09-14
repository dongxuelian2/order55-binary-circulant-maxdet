# Failed, timed-out, and withdrawn routes

The following outcomes are preserved so that collaborators do not repeat them
or mistake them for proofs.

- Entry-only ascent (200 and 4,000 basins) and hybrid row/column ascent (500
  basins) found no improvement.  These are negative heuristics.
- The 3^14 circulant scan ranked candidates with floating-point FFT scores and
  exactly checked 128 finalists.  It is not an exact proof of the unrestricted
  circulant optimum.
- Early integer and Boolean friendship SAT encodings returned `unknown` on
  timeouts.  The exact row-color congruence now excludes the all-norm-3 F7
  target independently.
- The corrected contingency-orbit replay has one certified UNSAT orbit and
  four timeout-unknown ten-row orbits.  Withdrawn pre-correction orbit claims
  must not be cited.
- Two OR-Tools downloads failed hash verification; CP-SAT was not installed
  or run.  No conclusion depends on it.
- A native Kissat timer interruption crashed the backend; this is an
  implementation failure, not UNSAT.
- A fresh triangle-component branch initially used an empty all-unit layer and
  returned zero.  It was corrected to the weighted all-unit layer and rerun
  with exact maxima.  The Q=126 K₄ forest norm-36 outside factor was likewise
  corrected from 198 to 189.  The corrected code and tests are the archived
  versions.

Any new timeout should be reported with solver, timeout, model revision, and
stdout; never silently convert `unknown` to `unsat`.
