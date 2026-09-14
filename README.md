# Binary circulant maximal determinant at order 55

This repository accompanies the manuscript **“Cyclotomic–Correlation Reduction for Binary Circulants: The Exact Order-55 Maximum.”** It contains the exact finite-enumeration code and the machine-readable proof data used in the order-55 classification.

## Main result

For binary circulant matrices of order 55,

```text
D_01(55) = 134694094094758395331307111329132
```

Every maximizing word has weight 28. The maximizers form a single affine-equivalence class under `i -> u i + t (mod 55)`, with 2,200 words and trivial stabilizer. A lexicographically least representative is

```text
0000000100011111011011001101011011010100011110101000111
```

The exact cyclotomic factors of this representative are

```text
N_5  = 131
N_11 = 434039
N_55 = 84603917386870862636041
```

and `28 * N_5 * N_11 * N_55` equals the displayed maximum.

## What the proof does

The proof is a finite computer-assisted proof with explicit mathematical coverage arguments. It does not rely on a heuristic search having explored “enough” words, and it does not require a second independently written lift program as a logical assumption.

The reduction is:

```text
all 2^55 binary words
  -> complement + exact weight bounds
  -> 5- and 11-fold cyclotomic constraints
  -> complete correlation-multiset recursion
  -> exact ordered correlation-profile screening
  -> 2,316 correlation targets
  -> 16,084 margin/correlation lift tasks
  -> exhaustive 5+6 meet-in-the-middle binary fibers
  -> one maximizing affine class
```

The manuscript proves the completeness of the fold enumeration, the correlation-multiset recursion, the target reduction, and the meet-in-the-middle lift. The code executes those finite reductions and records their outputs.

## Paper

The current single-file LaTeX manuscript is:

- [`paper/main.tex`](paper/main.tex)

All references are embedded in that file, so there is no separate BibTeX database. Compile it with two ordinary `pdflatex` passes:

```bash
cd paper
pdflatex main.tex
pdflatex main.tex
```

A detailed map from mathematical proof obligations to source files and certificate artifacts is in [`PROOF_MAP.md`](PROOF_MAP.md).

## Reproduce the stored certificate

Requirements:

- Python 3.11+
- SymPy
- Clang/LLVM (or another compatible compiler) with C++20 and unsigned `__int128` support

Install the package and run the verifier:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python scripts/verify_order55_global.py
python -m pytest
```

To regenerate the full finite enumeration before verification:

```bash
python scripts/build_order55_global.py --workers 8
python scripts/verify_order55_global.py
```

The optional `--full-lifts` mode repeats every binary lift using the complementary 6+5 split. It is a consistency check; the completeness of the 5+6 enumeration follows from the meet-in-the-middle lemma in the manuscript.

## Certificate summary

The committed `certificates/order55_global/global.json` records:

```text
status               COMPLETE
active weights        25, 26, 27
correlation targets   2,316
lift tasks             16,084
joined words           15,487,882,832
affine classes         1
maximizing words       2,200
```

The two primes used for exact signed CRT reconstruction are

```text
P_1 = 2305843009213696591
P_2 = 2305843009213697141
```

The SHA-256 digest of the certificate manifest recorded by `global.json` is

```text
b73a12c5b8502ac5cba20f5164ff351ba5df1a208d7cd2e67f45c359d7eedfc5
```

### Integrity note

The existing proof-critical generator/verifier sources and the stored certificate are intentionally kept byte-for-byte unchanged in this manuscript-only repository revision, because `source_hashes.json`, `hash_manifest.json`, and `global.json` bind those exact source bytes to the committed proof data. Historical internal names such as `*_audit.json` and the word `audit` in the verifier are therefore artifact nomenclature, not a statement that publication depends on an unfinished external audit. Changing those files would require regenerating and re-hashing the certificate.

`build.json` also contains the original machine’s absolute Windows paths. Those paths are nonsemantic provenance fields; the verifier resolves the committed source files by basename and checks their SHA-256 values. Compiler flags, source hashes, exact output files, counts, words, and determinant values are the relevant reproducibility data.

## Repository layout

- `paper/main.tex` — current journal-oriented manuscript, including references and the full weight table.
- `PROOF_MAP.md` — proof obligation → algorithm → artifact map.
- `scripts/build_order55_global.py` — regenerates the finite enumeration.
- `scripts/verify_order55_global.py` — verifies the stored enumeration, hashes, bounds, targets, witnesses, determinant and affine classification.
- `native/order55_profiles.cpp` — exhaustive 5- and 11-fold profile enumeration.
- `native/order55_correlation_screen.cpp` — exact ordered correlation-profile screening and CRT product reconstruction.
- `native/order55_torus_lift.cpp` — exhaustive meet-in-the-middle binary lift.
- `src/maxdet/order55.py` — exact circulant identities, determinant reconstruction and moment bounds.
- `src/maxdet/boxed_moment.py` — moment product bound with a spectral cap.
- `certificates/order55_global/` — immutable machine-readable proof data.
- `tests/` — focused exact-arithmetic tests.

## Order-15 third-root research track

The repository also contains the separate `order15_mu3/` research track for the order-15 third-root maximal-determinant problem.  This work is intentionally kept distinct from the proved order-55 binary-circulant theorem above.  The current rigorous order-15 record is accompanied by exact branch certificates through Q=150 and shellwise tail certificates for the only genuine remaining high-energy shells `Q=153,159,162,168`; `Q>=171` is below the record by trace stability.  See `order15_mu3/reports/CURRENT_STATUS.md` and `archive/stages/Q150-Q153.md` for the current frontier.

## Scope

The order-55 theorem concerns **55×55 binary circulant matrices**. It is not a claim about unrestricted binary matrices or unrestricted `{±1}` matrices. The repository separates exploratory discovery from the finite proof: only the exact threshold, mathematical reductions, exhaustive enumerations and exact determinant checks enter the theorem. The order-15 third-root material is a separate ongoing research program and is not part of that theorem.
