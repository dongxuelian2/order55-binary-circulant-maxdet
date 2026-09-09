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
