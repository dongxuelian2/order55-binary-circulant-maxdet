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

## Hadamard-ratio terminology

The exact rational ratios in the checkpoint are squared-determinant ratios:
`D_record/15^15 = 19365101568/30517578125` and
`U/15^15 = 60153061831/91552734375`.  The corresponding determinant-magnitude
ratios are display-only metadata,
`sqrt(D_record/15^15) ≈ 0.7965900126038639` and
`sqrt(U/15^15) ≈ 0.8105750078551661`, i.e. approximately 79.6590%–81.0575%
of the Hadamard magnitude bound.  These decimal square roots are never used as
proof inputs.

## What is proved and replayable

The proof spine is executable and integer/rational throughout:

1. `src/maxdet/mu3.py` implements Eisenstein arithmetic, exact determinants,
   dephasing, and the complete order-15 inner-product catalogue.
2. `verify_benchmark.py` recovers the record and checks the Gram determinant.
3. `verify_trace_stability.py` gives the strict energy restriction
   `Q=sum_{i<j}|G_ij|²≤168` for any strict counterexample.
4. The row-color congruence, the `F7` exclusion, the ternary-code independence
   bound, and the row/column intertwining identity reduce the search to four
   color partitions: `(15,0,0)`, `(14,1,0)`, `(13,2,0)`, `(13,1,1)`.  The same
   congruence gives `Q≡0 (mod 9)` for `(15,0,0)` and `(13,1,1)`, and
   `Q≡6 (mod 9)` for `(14,1,0)` and `(13,2,0)`; hence `Q≡3 (mod 9)` cannot
   occur.
5. Exact support/Schur certificates close the slices through Q=150.  The
   scripts named in `docs/CERTIFICATE_MAP.md` are the authoritative entry
   points; matching tests are under `tests/test_mu3_*.py`.
6. For Q≥153, stationary trace optimization is monotone on the genuine
   shells Q=153,159,162,168, and Q≥171 is already below the record.  The
   aggregate verifier also evaluates Q=156 and Q=165 conservatively; this does
   not change the Q=153 maximum or the certified integer upper.  The exact
   integer sieve uses `3^14` divisibility and the rational Eisenstein-norm test.

## What is not proved

There is no construction above `D_record`, no exact enumeration of the
high-energy Q=153,159,162,168 tail, and no independent full-fibre classification of
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
- `docs/ENVIRONMENT.md` — queried replay versions and clean installation.
- `order15_mu3/scripts/verify_aggregate_upper.py` — final aggregate sieve.
