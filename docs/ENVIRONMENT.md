# Environment record

This file records the environment used to validate the baseline repository.
It should be updated only when the reproducible experiment environment is
intentionally changed.

## Host

- OS: Windows 11 (10.0.26200)
- CPU: AMD Ryzen 7 7840H with Radeon 780M Graphics; 8 cores / 16 logical processors
- Python: 3.13.1
- Python executable: `C:\Users\29848\AppData\Local\Programs\Python\Python313\python.exe`

## Packages

The project installs the minimum requested runtime and test packages:

- NumPy: 2.5.3
- SciPy: 1.18.1
- SymPy: 1.14.0
- pytest: 9.1.1
- maximal-determinant: 0.1.0 (editable install)
- Numba: not installed; optional and intentionally skipped for this lightweight baseline

## Verification

- `pytest`: 22 passed in 1.90s
- `pip check`: `No broken requirements found.`
- Hadamard orders 1, 2, 4, 8: all exact-determinant tests passed
- Tiny random baseline: `n=8`, `samples=32`, `seed=0`; best fast score
  `255.99999999999994`, exact determinant `-256`
- Tiny normalized enumeration: known maxima for `n=1,2,3,4` all passed

The machine currently exposes a working Python 3.13 installation through the
`Python313` executable above. The generic `python`/`python3` commands resolve to
an unavailable Cygwin shim on this host, so commands in this project use the
virtual-environment interpreter explicitly.

## Order-55 global proof audit (2026-09-09)

- CPU SIMD probe: AVX2, BMI2, POPCNT, AVX512F, and AVX512VPOPCNTDQ available.
- Final native build: Clang C++20, `-O3 -march=native -Wall -Wextra`.
- Final suite: 30 passed; `pip check` reports no broken requirements.
- Original fixed-weight verifier: all seven strata pass, 29,332,216 covered
  zero-containing representatives in the original certificate.
- Production verifier: all mathematical bounds, hashes, exact norms, complete
  profile generators, listed task coverage, and all maximizing words PASS.
- Split-dependent 6+5 replay: all 16,084 listed margin tasks pass, with
  identical 15,487,882,832 joined-word counts and two normalized solutions;
  this is not an independent full fiber-union proof.
- Paper: 8-page PDF generated with local TeX Live 2024 and visually checked
  after Poppler rendering; no overfull-box or LaTeX warning remains.

## Order-15 μ₃ archival replay (2026-09-14)

The successful archived replay used `E:\maximal determinant\.venv\Scripts\python.exe`
with Python 3.13.1 and these installed packages:

- SymPy 1.14.0
- NumPy 2.5.3
- SciPy 1.18.1
- NetworkX 3.6.1
- z3-solver 5.1.0.0
- python-sat 1.9.dev15
- pytest 9.1.1

The full suite result for the handoff was 110 passed with 3 harmless NetworkX
warnings.  These versions were queried from the active virtual environment;
they are recorded here as provenance, not as mathematical inputs.

### Clean-environment installation

From a fresh checkout on Windows PowerShell:

```powershell
py -3.13 -m venv .venv
& .\.venv\Scripts\python.exe -m pip install -e ".[sat,graph]"
& .\.venv\Scripts\python.exe -m pytest -q
```

Then run `order15_mu3\scripts\replay_checkpoint.ps1` for the deterministic
certificate replay.  The generic `python` shim on the archival host was not
used; the explicit virtual-environment interpreter was.
