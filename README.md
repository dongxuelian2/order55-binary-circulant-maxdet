# Maximal determinant: exact binary circulant order 55

## Result

```text
D_01(55) = 134694094094758395331307111329132
```

All maximizing 55-bit words have weight 28 and form exactly one affine
class of size 2,200, with trivial stabilizer. Its lexicographically least
representative is

```text
0000000100011111011011001101011011010100011110101000111
```

This is a global theorem over all `2^55` binary circulants. The certificate
covers every weight pair, including equality, and classifies every maximizer.

- Paper: [Markdown](paper/order55_global/manuscript.md),
  [LaTeX](paper/order55_global/manuscript.tex),
  [PDF](output/pdf/order55_global.pdf).
- [Global certificate](certificates/order55_global/global.json) and
  [winner classification](certificates/order55_global/winner.json).
- [Independent verifier](scripts/verify_order55_global.py).
- [Source and novelty record](references/order55/SOURCES.md).
- [OEIS update draft](paper/order55_global/oeis_update_draft.md), not submitted.

## Reproduce

Install Python 3.11+ with this project (`python -m pip install -e .`) and
Clang with C++20 and unsigned `__int128` support. No external proof repository
or stochastic search is needed for the final certificate.

```powershell
& ./.venv/Scripts/python.exe scripts/build_order55_global.py
& ./.venv/Scripts/python.exe scripts/verify_order55_global.py --full-lifts
& ./.venv/Scripts/python.exe -m pytest
& ./.venv/Scripts/python.exe -m pip check
```

On other platforms use the installed Python interpreter instead of the
Windows venv path. The builder compiles native executables itself. For a
separate reproduction directory, pass `--output PATH`. The verifier accepts
that directory as its positional argument.

The ordinary verifier independently recomputes all bounds and norms,
recompiles the hash-bound generators, regenerates all folded and correlation
profiles, and checks the full task coverage. `--full-lifts` additionally
repeats every binary lift with a different 6+5 split. The original full
5+6 lift and this full independent split replay both completed successfully.

The proof reduces 38,629,684 formally allowed ordered correlation profiles
to 6,845,230 explicit profile evaluations, 2,316 correlation targets, and
16,084 complete margin tasks. The lift tests 15,487,882,832 joined words.
Two normalized solutions give the same affine class. The winner has three
independent exact checks: Fourier/CRT, SymPy integer determinant, and
cyclotomic resultant product.

The certificate stores exact source and output hashes. Preserve the bytes
of certificate files; their `.gitattributes` rules disable line-ending
conversion. Elapsed times and host-specific command paths may differ during
regeneration; mathematical counts, targets, and classification must agree.

## Supporting results and research history

The previously certified `D_circ(55,k)` values for `1 <= k <= 7` remain in
`certificates/order55_fixed_weight_1_7.json` and are discussed in the paper's
appendix. Run `scripts/verify_structured.py` to audit them.

Exploratory `scripts/order55_*.py` programs use scratch files under `results/`
and sometimes the separately retrieved ALETHEIA source. They are not
needed for the final proof. The final entry points are
`build_order55_global.py` and `verify_order55_global.py`.

The source checks on 2026-09-09 located no prior published exact order-55
determination in the sources inspected; this is not an absolute-priority
claim. No remote publication or OEIS submission has been performed.
