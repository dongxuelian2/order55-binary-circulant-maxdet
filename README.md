# Maximal Determinant Production Research

This repository contains exactness-aware code and certificates for determinant
searches over binary and ±1 matrices.

## Result

The production run of 2026-09-09 establishes the exact fixed-weight binary
circulant maxima `D_circ(55,k)` for `1 <= k <= 7`. The largest new structured
value is

```text
D_circ(55,7) = 570999857161750276477.
```

This is a structured exact theorem, not an unrestricted 55-by-55 binary
circulant result and not an unrestricted ±1 maximal determinant record. See
`paper/order55_fixed_weight_note.md` for the claim and proof, and
`certificates/order55_fixed_weight_1_7.json` for the machine-readable data.

## Reproduce

```powershell
clang++ -O3 -std=c++17 -Wall -Wextra -Wpedantic -pthread native/circulant_fixed_weight_exact.cpp -o native/circulant_fixed_weight_exact.exe
1..7 | ForEach-Object { native/circulant_fixed_weight_exact.exe --weight $_ --threads 15 }
& .\.venv\Scripts\python.exe scripts\verify_structured.py
& .\.venv\Scripts\python.exe -m pytest
& .\.venv\Scripts\python.exe -m pip check
```

The exact enumerator uses two finite fields and signed CRT for every candidate.
The independent verifier reconstructs the seven winning matrices and computes
their determinants with SymPy. Floating point is used only by the separate
bordered two-circulant heuristic searcher.

## Layout

- `src/maxdet/`: validated matrix, construction, and heuristic scoring primitives.
- `native/`: high-performance heuristic and exact exhaustive searchers.
- `scripts/verify_structured.py`: independent certificate verifier.
- `certificates/`: machine-readable exact results.
- `results/`: production-run discovery log and negative unrestricted evidence.
- `docs/target_selection.md`: literature-driven target selection.
- `references/SOURCES.md`: frozen upstream revisions and access dates.
- `paper/order55_fixed_weight_note.md`: technical note.

The original exact determinant infrastructure and Sylvester Hadamard fixtures
remain covered by the test suite.