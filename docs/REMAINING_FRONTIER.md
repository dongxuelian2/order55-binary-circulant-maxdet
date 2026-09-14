# Remaining mathematical frontier

The low and middle energy slices are covered by exact support, Schur, SAT/SMT,
or rational spectral certificates through Q=150.  The remaining gap is not a
known missing file or a failed replay: it is the sharpness of the continuous
high-energy bound.

## Tail shells

`Q=153,156,159,162,165,168` are the admissible shells in the trace tail.  The
stationary two-level product is largest at Q=153; the other shells are
strictly smaller.  `Q≥171` is below the published record by the independent
trace-stability certificate.  The aggregate verifier combines this tail with
the exact Q=144 and Q=150 branch closures, then applies `3^14` divisibility and
the rational Eisenstein norm test.

## What a successful continuation must do

1. Prove a structural restriction on Gram supports/phases at Q=153–168 that
   lowers the trace envelope, or enumerate every realizable support at those
   shells with exact arithmetic.
2. Keep both row and column μ₃ constraints; a one-sided Gram catalogue is not
   sufficient because the intertwining identity is a central obstruction.
3. Emit exact rational bounds and a deterministic replay script.  Floating
   ranking may be used to discover candidates, but not as the final claim.
4. Re-run `verify_aggregate_upper.py` and update the checkpoint only after the
   new bound is strict and its Eisenstein-norm sieve passes.

The best starting code is `verify_pre_q144_audit.py` (for component-catalogue
patterns), `verify_q150_boundary.py` (for Schur allocation patterns), and
`verify_trace_stability.py` (for the tail objective).
