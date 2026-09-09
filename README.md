# ±1 Maximal Determinant Baseline

This repository contains lightweight, reproducible infrastructure for
experiments on the maximal absolute determinant of `n × n` matrices with
entries in `{−1, +1}`.

The current stage is deliberately only baseline infrastructure. It provides
validated matrix primitives, exact determinant certification, a Sylvester
Hadamard fixture, tiny sanity checks, and a deterministic random-scoring
scaffold. It does not claim any new bound or start a research search.

`fast_determinant` and `fast_logabs_determinant` use NumPy floating-point
arithmetic for heuristic scoring only. A floating-point determinant is never a
mathematical certificate; claims and recorded candidates must use
`exact_determinant`, which uses SymPy exact arithmetic.

## Layout

- `src/maxdet/`: core primitives and the random baseline API.
- `tests/`: short pytest regression and sanity tests.
- `scripts/random_baseline.py`: command-line smoke baseline.
- `docs/environment.md`: the frozen local environment record.
- `data/` and `results/`: reserved for future small inputs and generated output.

## Setup and verification

On Windows, use the project interpreter explicitly:

```powershell
& .\.venv\Scripts\python.exe -m pip install -e .
& .\.venv\Scripts\python.exe -m pytest
& .\.venv\Scripts\python.exe -m pip check
& .\.venv\Scripts\python.exe scripts\random_baseline.py --n 8 --samples 1000 --seed 0
```

The random baseline accepts an optional `--output path.json`; JSON output
includes a UTC timestamp, Git commit, Python and package versions, parameters,
seed, fast score, exact determinant, and matrix.

Future work may add local search, switching, Gram-matrix methods, and
equivalence reduction. Those components are intentionally not implemented in
this baseline.
