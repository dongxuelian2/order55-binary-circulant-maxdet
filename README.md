# Maximum determinant of a 55×55 binary circulant matrix

This repository contains the exact finite-enumeration code and generated proof data accompanying the preprint **“The Maximum Determinant of a 55×55 Binary Circulant Matrix.”**

## Result

For binary circulant matrices of order 55,

```text
D_01(55) = 134694094094758395331307111329132
```

Every maximizing word has weight 28. The maximizers form a single affine-equivalence class under `i -> u i + t (mod 55)`, with 2,200 words and trivial stabilizer. A lexicographically least representative is

```text
0000000100011111011011001101011011010100011110101000111
```

The proof combines complement duality, exact spectral-moment bounds, the 5×11 cyclotomic decomposition, finite correlation-profile enumeration, and an exhaustive meet-in-the-middle lift of the remaining binary fibers. The mathematical coverage of the finite enumeration is proved in the paper; the code here executes those reductions and reproduces the reported counts and witnesses.

## Reproduce

Requirements:

- Python 3.11+
- Clang with C++20 and unsigned `__int128` support

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python scripts/verify_order55_global.py
```

To regenerate the finite enumeration before verification:

```bash
python scripts/build_order55_global.py --workers 8
python scripts/verify_order55_global.py
```

The optional `--full-lifts` verifier mode repeats the lift with the complementary 6+5 split as a consistency check. It is not required for the completeness proof.

## Repository layout

- `scripts/build_order55_global.py` — reconstructs the finite enumeration.
- `scripts/verify_order55_global.py` — checks the stored enumeration and final classification.
- `native/order55_profiles.cpp` — folded profile enumeration.
- `native/order55_correlation_screen.cpp` — exact correlation-profile screen.
- `native/order55_torus_lift.cpp` — meet-in-the-middle binary lift.
- `src/maxdet/order55.py` — exact circulant identities, determinant reconstruction, and moment bounds.
- `src/maxdet/boxed_moment.py` — product bound with a spectral cap.
- `certificates/order55_global/` — generated data used by the verifier.
- `tests/` — focused unit tests for the exact arithmetic routines.

The main branch is intentionally kept as a compact reproducibility repository. Historical audit material, exploratory searches, obsolete fixed-weight side projects, rendered paper build products, review-response drafts, and duplicate manuscript formats are omitted because they are not part of the proof of the order-55 theorem.
