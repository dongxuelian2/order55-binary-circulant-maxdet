# Exact maximal determinant over the third roots at order 15

This directory isolates the order-15 research program from the pre-existing
binary-circulant project in the repository root.

The matrix alphabet is represented logarithmically: `0,1,2` stand for
`1, omega, omega^2`, where `omega^2 + omega + 1 = 0`.  Exact certificate
arithmetic lives in `src/maxdet/mu3.py`; no floating-point result is used as a
mathematical claim.

Run the current benchmark verifier with:

```powershell
& ../.venv/Scripts/python.exe scripts/verify_benchmark.py
& ../.venv/Scripts/python.exe scripts/verify_bounds.py
```

See `reports/CURRENT_STATUS.md` for the concise interruption-safe state and
`docs/research_log.md` for chronology.

The archival handoff is maintained at repository level: see
`../docs/HANDOFF.md`, `../docs/COLLABORATOR_QUICKSTART.md`, and
`../results/handoff_checkpoint.json`.  The replay wrappers in this directory's
`scripts/replay_*.ps1` files keep the fast, checkpoint, and full commands
explicit.  The exact replay is deterministic but compute-bound; do not replace
it with a sampled or floating-point run when checking a certificate.
