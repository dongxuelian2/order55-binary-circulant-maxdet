# Response to the global-completeness review

Status: **MAJOR REVISION REQUIRED**.  The exact order-55 theorem claim and the
incumbent are retained.  The package is a review draft, not an arXiv-ready
submission, because the full order-55 fiber-union certificate has not yet been
independently reconstructed.

The production result under review is

```text
D_01(55) = 134694094094758395331307111329132
```

with one affine class of 2,200 weight-28 maximizers.  The canonical word,
its weight, complement, and complement determinant are unchanged:

```text
0000000100011111011011001101011011010100011110101000111
weight 28
complement 1111111011100000100100110010100100101011100001010111000
complement determinant 129883590734231309783760428781663
```

## Concern 1: the correlation screen was not independently checked

Agreed that replaying the production correlation executable would not be an
independent check.  We added a clean-room C++ evaluator under
`audit/completeness/cleanroom_profile_checker_v3.cpp`.  It reconstructs the
ordered profile permutations, fold tests, unit canonicalization, and the
finite-field profile product without importing the production verifier or
calling a production executable.

The resulting report is
`certificates/order55_global_v2/cleanroom_profile_report_v3.json`.  It gives:

- exact ordered-profile counts 407,277, 32,178,654, and 6,043,753 for
  weights 25, 26, and 27;
- exact canonical-profile counts 6,358, 822,567, and 155,980;
- 2,316 clean-room records at or above the incumbent;
- zero missing and zero extra records when compared as the set of
  `(weight, ordered correlation profile, exact product)` triples against the
  production target set.

The report prints `INDEPENDENT CORRELATION PROFILE CHECK: PASS`.  The clean
pruning table is stored as `pruned_profiles_v3.jsonl`; each row includes a
profile id, profile data, bound type, exact upper bound, incumbent, and the
pruned flag.

## Concern 2: the binary lift was replayed with the same implementation

Agreed.  The existing five-plus-six production run and its six-plus-five
replay are now described as a **split-dependent cross-check**.  They are not
called independent implementations in the revised manuscript.

The clean-room Python checker
`audit/completeness/independent_lift_checker.py` independently derives the
task descriptors from the target records and verifies:

- 2,316 targets expand to exactly 16,084 task descriptors;
- all 16,084 task descriptors occur exactly once in the eight task files;
- all 16,084 audit rows match the task partition;
- the two stored witnesses independently realize their task correlations and
  exact complementary determinants.

This establishes `listed_task_partition: PASS`, but it does **not** establish
the union of every full margin fiber.  Its report therefore records
`independent_full_fiber_union: NOT_COMPLETED`.  Closing this item requires a
second lift generator or a branch-level partition/union certificate that does
not reuse the same native lift implementation.  This is the outstanding
release blocker.

## Concern 3: `covered_words = 2**55` was a hard-coded assertion

Agreed.  The canonical `scripts/verify_order55_global.py` has been replaced by
the tested version whose final line reports

```text
domain_size = 1 << 55 (mathematical domain cardinality; not a derived coverage count)
```

The misleading `covered_words` label and `2**55` expression are gone from the
publication verifier and manuscript.  The numeric domain cardinality appears
only in the independent audit JSON as an explicitly labelled mathematical
domain size; it is not used as an empirical coverage total.

## Concern 4: no independently checkable global-completeness certificate

Agreed, and the revised package does not claim otherwise.  The audit hierarchy
now separates the following results:

| Component | Status | Meaning |
|---|---:|---|
| Correlation profile universe and target set | PASS | Clean-room profile evaluator matches the target set exactly. |
| Exact profile pruning records | PASS | Below-incumbent profile products are machine-readable and independently regenerated. |
| Lift task partition | PASS | All 16,084 listed tasks and audit rows match independently derived descriptors. |
| Lift witness checks | PASS | The two recorded witnesses match their correlations and determinants. |
| Affine canonicalization | PASS | Independent orbit size 2,200 and stabilizer size 1; both determinant values match. |
| Small-order end-to-end regression | PASS | Full brute-force/reduction checks at n=15 and n=21. |
| Random forward coverage | PASS | Exactly 1,000,000 samples; every sample is rigorously pruned or explicitly covered; unknown=0. |
| Independent union of every order-55 fiber | NOT COMPLETED | Requires the second generator or a branch-level union proof described above. |

The n=15/n=21 result is a regression of the reduction architecture against
brute force, not a substitute for the missing order-55 proof.  Likewise, the
one-million-sample result is a forward-coverage control, not an exhaustive
certificate.

## Revised publication decision

Outcome B applies.  The review bundle is useful for checking the exact theorem
claim, all current certificate objects, and the precise remaining gap, but it
must not be uploaded as a completed arXiv proof package.  The theorem has not
been weakened to a candidate statement; the release label is withheld solely
because the independently checkable full-fiber union is still missing.

The frozen production theorem commit is
`e77a1aa6afea16fe59a35f9c0de93f385f9b8781`.  The audit and packaging changes
are made on the review branch and are listed in the final provenance and hash
manifest.
