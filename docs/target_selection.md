# Production target selection, 2026-09-09

All determinant values below use the convention requested for this project:
an `n × n` matrix over `{−1,+1}`, with `D` denoting the normalized integer
`|det| / 2^(n-1)`.  The public two-circulant records were reconstructed from
the frozen Caltech verifier rather than transcribed from rounded tables.

| order | mod 4 | current normalized D | source/date | construction | upper bound | search assessment |
|---:|---:|---:|---|---|---|---|
| 51 | 3 | 17776121037665193653653203125 | arXiv:2608.22518v2, 2026-08-31 | bordered two-circulant, type 2 | `51^(51/2)/2^50` | new 3.1% record; unsuitable for immediate re-attack |
| 103 | 3 | 4982879878430668496143734745139892539874762668459055333727438310268207104 | Caltech commit `1f4804f`, 2026-09-07 | bordered two-circulant, type 1 | `103^(103/2)/2^102` | highest unrestricted priority: latest gain over predecessor was tiny |
| 107 | 3 | 25405109779472820154713362533412847329846084693257600588972842045966418068198 | arXiv:2608.22518v2, 2026-08-31 | bordered two-circulant, type 1 | `107^(107/2)/2^106` | recent 0.44% record; retained as benchmark |
| 111 | 3 | 139781659519566648611004967987048891981864344843163654605529677292458895404455125 | arXiv:2608.22518v2, 2026-08-31 | bordered two-circulant, type 3 | `111^(111/2)/2^110` | second priority: 1.26% recent gain |
| 115 | 3 | 826077036505625311168501810913750519509871482996572659976510854621009815010641902068 | Caltech commit `1f4804f`, 2026-09-07 | bordered two-circulant, type 1 | `115^(115/2)/2^114` | third priority: newer than the paper's type-3 record |
| 119 | 3 | 5230841367259683341310925504649737070470595687463364420259918451181874057850196588023699 | arXiv:2608.22518v2, 2026-08-31 | bordered two-circulant, type 3 | `119^(119/2)/2^118` | control target after three ansatz ceilings |

The scan covered orders 23 through about 150 at the level needed for target
selection.  Exact unrestricted maxima are already known through order 22;
orders divisible by four were deprioritized; older tables were treated as
stale wherever the 2026 repositories supplied replacements.  The Caltech v2
records at 51, 107, 111, 115, and 119 were explicitly excluded as discoveries.

The structured fallback was selected before the unrestricted run: binary
circulants of order 55 at fixed row weight.  Brent and Yedidia published the
unrestricted exact sequence through order 53, and the frozen ALETHEIA project
proves the unrestricted order-54 value.  No published order-55 fixed-weight
table or any of the exact integers in the resulting certificate was located
in searches of arXiv, OEIS, GitHub, and the cited literature.  Fixed weight is
a standard, non-ad-hoc stratum, it has a direct spectral formulation, and its
row-norm bound makes exact CRT exhaustion practical.

Selection ranking was: novelty, probability of improvement, computational
tractability, then certificate simplicity.  The executed targets and switching
evidence are recorded in `results/discovery_order55_fixed_weight.md`.
