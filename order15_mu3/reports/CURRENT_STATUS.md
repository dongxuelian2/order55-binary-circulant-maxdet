# Current status (archived 2026-09-14)

- **Current exact record:** recovered and exactly verified published `M15`,
  `det=604661760+241864704 omega` and
  `|det|^2 = 2^22 3^20 19 = 277868041444786176`.
- **Rigorous upper bound:** `287710229992756239`, from the full unmocked
  pre-Q144 replay, exact Q=144/Q=150 branch certificates, the Q≥153 trace
  envelope, `3^14` divisibility, and the Eisenstein norm sieve.
- **Exact maximality proved:** no.
- **Strongest verified theorem:** every strict counterexample has `Q<=168`
  after the trace-stability replay, and the four row/column color partitions
  are `(15,0,0)`, `(14,1,0)`, `(13,2,0)`, `(13,1,1)`.  Exact branch audits
  close all slices below the Q=153 tail envelope.
- **Important replay artifacts:** `scripts/verify_color_15_energy90.py`,
  `scripts/verify_color_15_energy99.py`,
  `scripts/verify_color_15_energy108.py`,
  `scripts/verify_q96_boundary.py`, `scripts/verify_q105_boundary.py`,
  `scripts/verify_q114_boundary.py`, and `scripts/verify_refined_upper.py`.
- **Remaining gap:** exact maximality is open; the upper/lower ratio is
  `60153061831/58095304704 ≈ 1.03542036895` and the additive squared gap is
  `9842188547970063`.  The open mathematical
  frontier is the unenumerated high-energy tail represented by the Q=153
  trace envelope.  The color congruence leaves genuine shells Q=153,159,162,168;
  the aggregate verifier also evaluates Q=156 and Q=165 conservatively.  Q>=171
  is below the record by trace stability.
- **Completed Q=117 partial searches:** every mixed isolated and mixed
  no-isolate partition is empty above the record; the no-isolate thirteen-unit
  layer has 10,329 exact profiles and no arithmetic survivor.  These are
  preserved as partial certificates, not promoted to a global enumeration.
