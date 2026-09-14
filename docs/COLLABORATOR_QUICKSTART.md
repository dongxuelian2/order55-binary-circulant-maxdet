# Collaborator quickstart

1. Read `results/handoff_checkpoint.json` and `docs/HANDOFF.md`.
2. Create a clean Python 3.11+ environment and install the project (plus
   `.[sat,graph]` if SAT or graph scripts are needed).
3. Run `order15_mu3/scripts/replay_fast.ps1`; confirm exact benchmark and
   representative Q=126 tests pass.
4. Run `order15_mu3/scripts/replay_checkpoint.ps1` on an uninterrupted host.
   Save stdout and compare the returned exact integers with the checkpoint.
5. Before changing mathematics, read `docs/PROOF_DEPENDENCY_GRAPH.md` and
   `docs/FAILED_AND_TIMEOUT_ROUTES.md`.  Preserve the four color partitions and
   both-sided row/column constraints.

The single best file to open first is
`docs/HANDOFF.md`; the single best executable entry point is
`order15_mu3/scripts/verify_aggregate_upper.py`.
