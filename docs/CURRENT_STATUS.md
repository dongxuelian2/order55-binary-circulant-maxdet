# Current status — order 15 over μ₃

As of 2026-09-14:

| Item | Status |
|---|---|
| Published lower bound | `277868041444786176 = 2^22·3^20·19` |
| Exact benchmark replay | Passed (`verify_benchmark.py`, exit 0) |
| Full unmocked pre-Q144 replay | Passed (`verify_pre_q144_audit.py`, exit 0) |
| Q=144/Q=150 branch certificates | Implemented and tested |
| Q≥153 trace envelope | Exact rational replay; Q=153 is the largest checked tail shell |
| Current integer upper | `287710229992756239` after `3^14` + Eisenstein norm sieve (verify checkpoint) |
| Exact maximality | **Open** |
| Best next target | Replace the Q=153–168 continuous envelope by a realizability or phase certificate |

The upper/lower ratio, Hadamard ratios, exact rational envelope, and replay
metadata are intentionally duplicated in
`results/handoff_checkpoint.json`; update this page only when that file and the
verifier output change together.
