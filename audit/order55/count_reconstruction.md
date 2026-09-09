# Count reconstruction

All counts below were reconstructed from the committed certificate and then
checked by audit-only code.  A count is accepted only when its definition and
its independent recomputation agree.

## Weight and fold counts

| Object | k=25 | k=26 | k=27 |
|---|---:|---:|---:|
| Maximum allowed half-square sum | 3338 | 3919 | 4570 |
| m=5 folded margin survivors | 17 | 27 | 26 |
| m=11 folded margin survivors | 321 | 560 | 660 |

The weight screen has exactly 27 rows.  Every row has the expected complement
weight, upper-floor formula, and exclusion flag; only 25, 26, and 27 remain
eligible for the global maximum.

## Correlation profile partitions

The independent profile screen uses a different traversal architecture from
the production recursive enumerator.  It starts each sorted multiset and
visits every distinct permutation with `std::next_permutation`.

| k | part | ordered | folded pass | canonical exact | >= maximum | targets |
|---:|---:|---:|---:|---:|---:|---:|
| 25 | 2 | 2,925 | 1,065 | 55 | 55 | 55 |
| 26 | 3 | 807,300 | 440,990 | 22,094 | 0 | 0 |
| 26 | 5 | 8,775 | 5,440 | 279 | 279 | 279 |
| 26 | 6 | 27 | 25 | 2 | 2 | 2 |
| 27 | 1 | 5,920,200 | 3,051,280 | 152,721 | 0 | 0 |
| 27 | 3 | 105,300 | 56,165 | 2,821 | 1,964 | 1,964 |
| 27 | 4 | 702 | 300 | 16 | 16 | 16 |
| 27 | 5 | 1 | 0 | 0 | 0 | 0 |
| **total** |  | **6,845,230** | **3,555,265** | **177,988** | **2,316** | **2,316** |

For partitions where the production interval screen stops before an exact CRT
product is evaluated, its diagnostic `max_absolute_profile_product` is zero.
The independent screen evaluates the exact product for every canonical profile;
the comparison therefore requires all traversal and target counts, and compares
the target sets exactly rather than comparing that interval-dependent diagnostic.

## Lift reconstruction

| part | tasks | joined words | solutions |
|---:|---:|---:|---:|
| 0 | 2,011 | 1,947,432,320 | 0 |
| 1 | 2,011 | 1,922,675,330 | 1 |
| 2 | 2,011 | 1,951,395,543 | 1 |
| 3 | 2,011 | 1,932,139,664 | 0 |
| 4 | 2,010 | 1,941,341,181 | 0 |
| 5 | 2,010 | 1,920,036,010 | 0 |
| 6 | 2,010 | 1,944,066,646 | 0 |
| 7 | 2,010 | 1,928,796,138 | 0 |
| **total** | **16,084** | **15,487,882,832** | **2** |

The clean-room split-6 solver reproduces every row of this table and the two
normalized output words.  Its audit output is hashed in
`independent_lift_audit.json`.

## Winner counts

The independent winner checker obtains:

* high word weight 28 and determinant
  `134694094094758395331307111329132`;
* complement weight 27 and determinant
  `129883590734231309783760428781663`;
* complement ratio on both sides
  `3636740540558476673945292005886564`;
* affine stabilizer 1 and orbit size 2,200;
* exactly one affine class in the committed winner certificate.

These counts are also checked by direct 55 by 55 Bareiss elimination, exact
cyclotomic resultants, modular Gaussian elimination, and Fourier products.
