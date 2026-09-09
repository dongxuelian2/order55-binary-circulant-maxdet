# Coverage audit

The audit treats every reduction as a map whose domain and image are explicit.
The final lift stage was recomputed, rather than accepted from a counter file.

```text
all 2^55 binary words
        |
        v  complement identity and determinant symmetry
weights 0..27 (zero endpoints removed)
        |
        v  first/fourth-moment plus rational spectral caps
active weights 25, 26, 27
        |
        v  complete capped 5- and 11-margin enumeration
folded margin survivors and folded-correlation signatures
        |
        v  sorted half-correlation multisets and unit canonicalization
exact correlation-profile targets (2316)
        |
        v  independently rebuilt row/column margin map
16,084 lift tasks
        |
        v  all compatible binary column choices, exact shifts 1..27
15,487,882,832 joined words -> 2 normalized low words
        |
        v  complement and affine orbit
one affine class, 2,200 maximizing words
```

## Domain accounting

| Layer | Required domain or count | Independent check |
|---|---:|---|
| Binary inputs | `2^55 = 36,028,797,016,963,968` | exact finite domain identity |
| Complement reduction | weights 0 through 27 | rederived Fourier complement identity |
| Weight rows | 27 rows, weights 1 through 27 | rational first/fourth-moment recomputation |
| Active weights | 25, 26, 27 | independent bound floors and thresholds |
| Fold profiles | m=5 and m=11 for each active weight | clean-room capped margin enumeration and exact resultants |
| Correlation profile universe | all stored square-bounded multisets | independent multiset recursion and multinomial counts |
| Exhaustive profile partitions | 8 partitions | clean-room `next_permutation` replay |
| Correlation targets | 2,316 | exact target-set equality in all 8 partitions |
| Lift tasks | 16,084 | independently rebuilt target-to-task topology and per-part counts |
| Lift joins | 15,487,882,832 | clean-room split-6 typed-key solver |
| Normalized solutions | 2 | independent output set equality and determinant checks |
| Affine orbit | 2,200 words | unit/translation orbit enumeration, stabilizer 1 |

## No hidden quotient assumption

The sorted correlation profile is only a quotient for screening.  It is not
treated as a binary realization.  Every canonical target is mapped back to all
compatible oriented 5- and 11-margin signatures.  The lift input records the
resulting task identifiers, row margins, column margins, and all 55 target
correlations.  The clean-room Python checker independently reconstructs this
orientation map and verifies each task's dimensions and identifiers.

The independent C++ lift solver then enumerates column masks for each task.  It
uses a seven-coordinate typed signature in an unordered map, joins every
matching left/right partial word, and tests every shift 1 through 27 after the
join.  Thus the 11 and 22 shift coordinates are not silently assumed to imply
the other correlations.

## Partition totals

The eight exhaustive profile partitions have ordered counts
`2,925 + 807,300 + 8,775 + 27 + 5,920,200 + 105,300 + 702 + 1`
= `6,845,230`.  Their independent replay reproduces the stored fold-pass,
canonical, and threshold counts and returns exactly 2,316 target records.

The eight lift parts have task counts
`2011, 2011, 2011, 2011, 2010, 2010, 2010, 2010`, summing to 16,084.  No
part is absent and no task is duplicated in the committed task list.

## Coverage conclusion

No unexplained gap remains in the active search.  The only non-exhaustive
operation is the one-million-word random forward test in the Python audit; it
is a diagnostic, not a proof step.  The proof-relevant profile and lift stages
were checked exhaustively over their declared finite domains.
