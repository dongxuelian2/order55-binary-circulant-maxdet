# Certificate map

| Claim or stage | Exact entry point | Matching test/status |
|---|---|---|
| μ₃ arithmetic, dephasing, catalogue | `src/maxdet/mu3.py` | `tests/test_mu3.py` |
| Published M15 lower bound | `order15_mu3/scripts/verify_benchmark.py` | benchmark test |
| Coarse and trace bounds | `verify_bounds.py`, `verify_trace_stability.py` | `test_mu3_trace.py` |
| Friendship and color reduction | `verify_friendship_bound.py`, `verify_color_partition_bounds.py` | color-bound tests |
| Q=45–81 all-same layers | `verify_color_15_low_energy.py`, `verify_color_15_energy81.py`, SAT scripts | `test_mu3_color_15*.py` |
| Q=60/69/78 mixed obstructions | `verify_color_14_1_twin_bound.py`, `verify_q78_color_orbit_obstruction.py` | corresponding tests |
| Q=87–117 boundary audits | `verify_q87_boundary.py`, `verify_q96_boundary.py`, `verify_q105_boundary.py`, `verify_q114_boundary.py`, `verify_q117_boundary.py` | boundary tests |
| Q=123 and Q=126 exact support branches | `verify_q123_boundary.py`, `verify_q126_boundary.py`, `verify_q126_*` | exact-certificate tests |
| Q=132/135/141 refinements | `verify_q132_boundary.py`, `verify_q135_boundary.py`, `verify_q141_boundary.py`, `verify_q141_energy36_minimal_cross.py` | matching tests |
| Full Q<144 aggregate replay | `verify_pre_q144_audit.py` | `test_mu3_pre_q144_audit.py` (long) |
| Q=144 and Q=150 closures | `verify_q144_boundary.py`, `verify_q150_boundary.py` | matching tests |
| Final integer/Eisenstein sieve | `verify_aggregate_upper.py` | `test_mu3_aggregate_upper.py` (long) |

The stage scripts return JSON to stdout.  No generated JSON is silently treated
as a certificate unless the corresponding script assertion succeeds.
