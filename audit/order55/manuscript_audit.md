# Manuscript and novelty audit

## Scope checked

The order-55 manuscript and its linked documentation were checked against the
claimed value, the high/low witness words, the one-class/2200-word statement,
the Fourier proof route, the weight and fourth-moment bounds, the 5/11 folds,
the correlation target count, and the lift coverage tables.

No mathematical contradiction was found.  The statement of the result is
consistent with the committed winner certificate, and the manuscript explains
why an affine orbit has size 2200 when the stabilizer is trivial.

## Required wording corrections

1. The manuscript calls the production lift replay a separate or independent
   verifier.  `scripts/verify_order55_global.py` recompiles and invokes the same
   native lift source used by the generator; its 6+5 run changes the split but
   not the implementation.  It should be described as a production replay
   under an alternate split, or the clean-room audit implementation should be
   cited as the independent verifier.
2. Reproducibility claims should state that the certificate binds seven source
   files and 71 certificate components, but not the compiler binary, Python
   environment, native executable, or all project configuration.  Add those
   digests in a revised release.
3. The historical incumbent-search commit should be labelled as historical
   provenance, distinct from the final frozen global revision.

These corrections are disclosure revisions.  They do not invalidate the
maximum or the uniqueness result because the missing independence has been
replaced here by clean-room exhaustive profile and lift replays.

## Reproducibility commands

The important audit commands are:

```text
.venv\\Scripts\\python.exe audit/order55/verify_winner_independent.py --random-count 100
.venv\\Scripts\\python.exe audit/order55/verify_global_independent.py --random-count 1000000
.venv\\Scripts\\python.exe audit/order55/run_independent_profiles.py --workers 8
.venv\\Scripts\\python.exe audit/order55/run_independent_lifts.py --workers 8
.venv\\Scripts\\python.exe audit/order55/small_order_regression.py
.venv\\Scripts\\python.exe -m pytest
.venv\\Scripts\\python.exe -m pip check
```

They produce PASS outputs and the JSON evidence files in this directory.

## Scoped novelty spot-check

A targeted web search for the exact 33-digit value, `D_01(55)`, the canonical
word, and binary-circulant order-55 maxima did not return an exact prior match.
It did return the known binary-circulant tables through order 25 in Richard
Brent's computation paper and the associated maxdet page:

* [Computation of Maximal Determinants of Binary Circulant Matrices](https://maths-people.anu.edu.au/~brent/pd/rpb270v1.pdf)
* [Richard Brent's Hadamard maximal determinant page](https://maths-people.anu.edu.au/~brent/maxdet/)

This is only a scoped novelty screen, not a literature proof.  The audit makes
no claim of priority.
