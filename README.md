# Maximum determinants of order-55 circulants over both binary alphabets

The exact global maxima are:

~~~text
D_01(55)      = 134694094094758395331307111329132
D_pm_circ(55) = 297532404423965431213849494795385191907933028352
S_pm(55)      = 16516366298178510386365752382353
~~~

The sign normalization is exactly D_pm_circ(55) = 2^54 S_pm(55).
The 0/1 maximum has 2,200 weight-28 maximizers in one affine orbit.
The sign maximum has 4,400 maximizers, of weights 23 and 32, in two
affine orbits of size 2,200 with trivial stabilizers. Global negation
pairs those orbits into one equivalence class.

Canonical lower-weight binary encoding of the sign winner:

~~~text
0000000000100110110100111001001110001010101100010010111
~~~

For x=2a-1 with binary weight k,
det C(x)=2^54 (2k-55) det C(a)/k. The sign optimum differs from the
binary optimum; every sign-dependent pruning threshold was regenerated.

## Proof and audit

The sign proof reduces 9,531,797,706 formally allowed ordered correlations
to 731,378,844 explicit profile evaluations, 17,835 targets, and 151,772
binary margin tasks. Production and independent full lift implementations
each inspect 113,997,335,580 joined words and agree task by task.

Independent arithmetic, multiset, FKM-necklace fold, profile, fiber-union,
and affine audits pass. The raw independent algebraic profile set has
255 additional indefinite spectra, each excluded by an explicit rational
negative-eigenvalue upper bound. The remaining target sets agree exactly.
Full sign regressions at 15 and 21 compare all words and all equality cases.

The original 0/1 certificates are preserved. The independent full binary
fiber supplement now resolves the outstanding condition in the historical
v2 review snapshot; all 16,084 original tasks agree with production.

- Joint paper: [PDF](output/pdf/order55_joint.pdf),
  [formal LaTeX](paper/order55_global/arxiv/main.tex),
  [Markdown](paper/order55_global/manuscript.md).
- Sign [global certificate](certificates/order55_sign/global.json),
  [winner](certificates/order55_sign/winner.json),
  [classification](certificates/order55_sign/classification.json).
- Sign [profile audit](certificates/order55_sign/cleanroom_profile_report.json)
  and [full fiber audit](certificates/order55_sign/cleanroom_lift_report.json).
- Original [binary certificate](certificates/order55_global/global.json).
- Completed [binary full-fiber supplement](certificates/order55_sign/binary_companion/full_fiber_report.json).
- [Source check](research/order55_sign/sources.md) and
  [final novelty record](certificates/order55_sign/novelty_check.json).

## Verify and reproduce

Use Python 3.11+ with SymPy and Clang/GCC supporting C++20 and unsigned
128-bit integers. The proof scripts select this checkout's source tree.

~~~powershell
python scripts/verify_order55_sign.py
python scripts/verify_order55_sign.py --replay
python scripts/render_order55_paper.py
~~~

The first command checks the frozen sources and assets, exact winner,
weight bounds, target sets and negative-spectrum discharges, full task
records, and all maximizing words. The replay command recompiles and
regenerates the complete sign pipeline, both independent audits, and
full small-order regression. It requires no heuristic search or network.
Paths and elapsed times can differ; mathematical counts and word sets
must agree.

The existing binary verifier and fixed-weight-1..7 certificates remain
available. The new sign global upper bound does not depend on the old
binary global theorem.

The work is isolated on branch sign-circulant-order55. No remote push,
publication, or OEIS submission was performed. The OEIS drafts for both
sign sequences accompany the paper; the A086432 draft is preserved.
