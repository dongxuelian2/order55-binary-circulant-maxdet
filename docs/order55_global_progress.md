# Order 55 global proof: completed

Starting commit: decd33b5ca8e428bf3dc2cc9f9db43938916c8cb.
Initial clean checkout and baseline: 26 tests passed.

The global result is 134694094094758395331307111329132. There is one
maximizing affine class, with trivial stabilizer and 2200 words, all weight
28. The final proof, complete classification, and independent verifier are
under `certificates/order55_global/`, `paper/order55_global/`, and
`scripts/verify_order55_global.py`.

The exploratory run found exact complementary incumbents
91089013395259690896562834081004,
95016526411094383114560322223123,
103755181499155107923236157710676, then the final value through exact MITM.
The third incumbent had no improving move in the executed complete
three-swap audit. This was only search evidence. Global proof came from
correlation and binary-lift coverage.

At the earlier threshold, weight 24 was exhausted through 296010 profiles
and eight margin fibers; weight 25 through 432505269 profiles and 45888
margin tasks, with 38147109272 joined words and zero solutions. Several
larger exploratory correlation jobs were stopped only after a new exact
capped-moment bound excluded their entire multiset. Their incomplete output
is scratch data and is not used by the final certificate.

The final higher threshold permits a much smaller, fresh complete proof:
6845230 ordered profiles explicitly enumerated, 2316 correlation targets,
16084 margin tasks, and 15487882832 joined words. That entire proof was
regenerated from a clean C++ build; the independent verifier recomputed its
bounds, exact norms, profile coverage, and all lifts using a different split.

The final proof does not rely on live computation, heuristic scores, old
partial job outputs, or external repositories. Scratch scripts remain for
research provenance; use the final builder and verifier for reproduction.
