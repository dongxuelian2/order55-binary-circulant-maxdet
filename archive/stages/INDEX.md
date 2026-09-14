# Stage index

| Stage | Archived status | Main evidence |
|---|---|---|
| Q=45–81 | Exact all-same low-energy obstructions; Q=81 reduced CNF UNSAT | `verify_color_15_low_energy.py`, `verify_color_15_energy81.py`, `sat_q81_reduced.py` |
| Q=90 | Exact top/second energy-9 obstructions and all-same cap | `verify_color_15_energy90.py`, `verify_q87_q90_*` |
| Q=99–108 | Exact all-same support layers; no arithmetic survivor above record | `verify_color_15_energy99.py`, `verify_color_15_energy108.py` |
| Q=117 | Mixed branches closed; thirteen-unit no-isolate layer enumerated; isolated layer partial | `verify_q117_boundary.py`, `verify_color_15_energy117.py` |
| Q=123 | Exact color/Schur allocations below Q=126 | `verify_q123_boundary.py` |
| Q=126 | Component/forest/K₄/K₅/triangle branches exact; all below Q=153 | `verify_q126_*` |
| Q=132–141 | Boundary and minimal-cross refinements | `verify_q132_boundary.py`, `verify_q135_boundary.py`, `verify_q141_*` |
| Q=144 | Exact boundary below Q=150 | `verify_q144_boundary.py` |
| Q=150–153 | Q=150 exact boundary below Q=153 tail; Q≥153 remains trace-only | `verify_q150_boundary.py`, `verify_aggregate_upper.py` |

Individual notes preserve the distinction between an exact closure, a
rational spectral cap, and an unresolved solver/search route.
