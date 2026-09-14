# Reproduction

All commands below run from the repository root on Windows PowerShell.  The
checked-in wrappers select `.venv\Scripts\python.exe` when present and fall
back to `python`.

## Fast smoke

```powershell
& .\order15_mu3\scripts\replay_fast.ps1
```

This exercises exact ring arithmetic, trace bounds, representative Q=126
certificates, the Q=141 residual certificate, and the Q=150 boundary.

## Checkpoint replay

```powershell
& .\order15_mu3\scripts\replay_checkpoint.ps1
```

This runs the benchmark, coarse bounds, the complete pre-Q144 audit, Q=144,
Q=150, and the aggregate sieve.  It is deterministic but compute-bound; allow
several minutes and do not add a timeout.

## Full test/replay path

```powershell
& .\order15_mu3\scripts\replay_full.ps1
```

The final `pytest` invocation currently collects 110 tests.  The full suite
includes optional SAT/SMT model-shape checks; solver `unknown_timeout` values
are data, not proofs.  A clean replay must report exit code 0 for each exact
certificate command.

## Direct Python entry points

```powershell
& .\.venv\Scripts\python.exe -m order15_mu3.scripts.verify_benchmark
& .\.venv\Scripts\python.exe -m order15_mu3.scripts.verify_pre_q144_audit
& .\.venv\Scripts\python.exe -m order15_mu3.scripts.verify_aggregate_upper
& .\.venv\Scripts\python.exe -m pytest -q
```

The pre-Q144 and aggregate commands print large JSON objects.  Preserve the
stdout when making a release record; the compact summary is in
`results/handoff_checkpoint.json`.
