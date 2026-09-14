# Order-15 μ₃ environment record

Snapshot: 2026-09-14, Windows PowerShell, repository `E:\maximal determinant`.

- Python interpreter used for exact replay: `E:\maximal determinant\.venv\Scripts\python.exe`.
- Project requires Python ≥3.11 (see `pyproject.toml`).
- Core dependencies: NumPy 2.5.3, SciPy 1.18.1, SymPy 1.14.0, pytest 9.1.1.
- Optional research dependencies present in the replay environment:
  `z3-solver 5.1.0.0`, `python-sat 1.9.dev15`, and NetworkX 3.6.1.
- Python version: 3.13.1.  The replay emits a harmless NetworkX graph-hashing
  warning under current versions; exact values are unaffected.
- No credentials, API keys, or external proof repositories are needed.
- `.venv`, caches, build products, and transient generated results are
  excluded from the handoff.  Only `results/handoff_checkpoint.json` is
  promoted under the root `results/` directory.

For a fresh environment, install the project and optional packages, then run
`order15_mu3/scripts/replay_fast.ps1` before attempting the long replay.
