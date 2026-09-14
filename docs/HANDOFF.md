# Order-15 μ₃ maximal-determinant handoff

Snapshot date: 2026-09-14 (Asia/Shanghai).  This is a reproducibility and
collaborator handoff for the exact order-15 third-root-of-unity investigation.
It does not claim a completed maximality proof.

## One-minute status

The published construction `M15` is reproduced exactly:

```text
det(M15) = 604661760 + 241864704 ω
D_record  = |det(M15)|² = 277868041444786176 = 2²²·3²⁰·19.
```

The complete unmocked `Q<144` replay (`verify_pre_q144_audit.py`) exited 0.
The aggregate verifier closes the exact Q=144 and Q=150 branches and applies
the two-level trace envelope to Q≥153.  Its exact integer/Eisenstein-norm
sieve gives the current rigorous upper
`287710229992756239` (nine non-norm multiples of `3^14` are skipped).  Thus
exact maximality remains open, with a squared gap of about 3.54%.

The direct literature audit found the record in Nuñez Ponasso's 2025 LAA
paper and February 2026 dissertation; both leave order 15 unproved.  See
`order15_mu3/docs/literature.md` for the bounded source audit.

## What is proved and replayable

The proof spine is executable and integer/rational throughout:

1. `src/maxdet/mu3.py` implements Eisenstein arithmetic, exact determinants,
   dephasing, and the complete order-15 inner-product catalogue.
2. `verify_benchmark.py` recovers the record and checks the Gram determinant.
3. `verify_trace_stability.py` gives the strict energy restriction
   `Q=sum_{i<j}|G_ij|²≤168` for any strict counterexample.
4. The row-color congruence, the `F7` exclusion, the ternary-code independence
   bound, and the row/column intertwining identity reduce the search to four
   color partitions: `(15,0,0)`, `(14,1,0)`, `(13,2,0)`, `(13,1,1)`.
5. Exact support/Schur certificates close the slices through Q=150.  The
   scripts named in `docs/CERTIFICATE_MAP.md` are the authoritative entry
   points; matching tests are under `tests/test_mu3_*.py`.
6. For Q≥153, stationary trace optimization is monotone on the checked
   shells Q=153,156,159,162,165,168, and Q≥171 is already below the record.
   The exact integer sieve uses `3^14` divisibility and the rational
   Eisenstein-norm test.

## What is not proved

There is no construction above `D_record`, no exact enumeration of the
high-energy Q=153–168 tail, and no independent full-fibre classification of
all matrices.  SAT/SMT timeouts are recorded as unknown rather than UNSAT.
The Q=117 isolated thirteen-unit layer is a historical partial search, not a
global certificate.  A collaborator must not promote any of these statements
to “maximal” without a new proof closing the tail or improving the bound.

## First files to read

- `results/handoff_checkpoint.json` — machine-readable constants and replay
  status.
- `docs/CURRENT_STATUS.md` — concise human status.
- `docs/REPRODUCTION.md` — commands and expected runtimes.
- `docs/CERTIFICATE_MAP.md` — theorem-to-script/test map.
- `docs/REMAINING_FRONTIER.md` — the best next mathematical attack.
- `order15_mu3/scripts/verify_aggregate_upper.py` — final aggregate sieve.
