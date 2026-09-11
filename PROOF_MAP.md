# Proof map for the order-55 theorem

This document connects each mathematical reduction in the manuscript to the exact program that executes it and to the committed machine-readable artifact that records the result. The repository is intended to make a computer-assisted proof inspectable: the manuscript supplies the completeness arguments; the programs execute the resulting finite enumerations.

## 1. Claimed theorem

For a binary word `a` of length 55, let `C(a)` be the corresponding circulant matrix. The claimed exact value is

```text
D_01(55) = 134694094094758395331307111329132.
```

Every maximizing word has weight 28, and all maximizing words form one affine orbit under `i -> u i + t (mod 55)`. The orbit has size 2,200 and trivial stabilizer.

A canonical representative is

```text
0000000100011111011011001101011011010100011110101000111
```

The two exact finite-field residues, the cyclotomic factors, the autocorrelation vector and the fold vectors are stored in `certificates/order55_global/winner.json`.

## 2. Logical structure

The global proof has the following implications.

```text
binary word with determinant >= M
        |
        | complement duality + exact moment bounds
        v
lower weight k in {25,26,27} with a short square-sum window
        |
        | prime-factor sector moments + fold norm thresholds
        v
realizable 5-fold and 11-fold correlation signatures
        |
        | complete sorted-multiset recursion + safe spectral cap
        v
ordered correlation profile surviving exact screening
        |
        | unit canonicalization
        v
one of 2,316 correlation targets
        |
        | all compatible fold orientations
        v
one of 16,084 binary lift tasks
        |
        | complete 5+6 meet-in-the-middle fiber enumeration
        v
two normalized lifted witnesses
        |
        | exact determinant + affine canonicalization
        v
one maximizing affine class of size 2,200
```

The first, third and fifth coverage arrows are stated and proved explicitly in `paper/main.tex`; the certificate records their finite outputs.

## 3. Weight screen

**Mathematical input.** For lower weight `k <= 27`, complement duality gives

```text
det C(1-a) = (55-k) * product(q_r),  r=1,...,27.
```

The first two power sums of the nonzero squared Fourier magnitudes `q_r` are fixed by the weight and by the square sum of the periodic autocorrelations. The product is bounded by the finite two-value extremum described in the manuscript.

**Implementation.** `src/maxdet/order55.py` implements the rational AM-GM/moment bounds. `scripts/verify_order55_global.py` recomputes every weight bound from exact rational arithmetic.

**Artifact.** `certificates/order55_global/weight_screen.json` contains all lower weights `1,...,27` and their integer upper bounds. It proves that `k <= 24` is below the incumbent. The only active weights are:

```text
k=25: sum c_s^2 in [3336,3338]
k=26: sum c_s^2 in [3913,3919]
k=27: sum c_s^2 in [4563,4570]
```

The full table is reproduced in the manuscript appendix.

## 4. Prime-factor folds

For `m=5` and `m=11`, the fold vector is

```text
x_j = sum a_t  over t == j (mod m).
```

Its autocorrelation is the corresponding residue-class sum of the full correlation vector. Parseval gives exact first and second moments of the order-`m` Fourier sector. AM-GM on the remaining sectors gives a necessary cyclotomic-norm threshold.

**Completeness argument.** Every fold vector has an affine representative whose minimum occurs in coordinate zero. `native/order55_profiles.cpp` recursively enumerates all bounded integer compositions of weight `k`, pruning only when a necessary remaining-sum or square-sum condition fails, and then canonicalizes against all affine maps on `Z/mZ`. The manuscript proves that this visits one representative of every feasible orbit.

**Exact norm.** The program evaluates the fold norm modulo

```text
P_1 = 2305843009213696591.
```

For `m <= 11` and `k <= 27`, the elementary bound `N_m <= k^(m-1) < P_1` makes the modular residue the unique ordinary integer norm.

**Artifacts.** `profiles_m*_k*.jsonl`, `foldcorr_m*_k*.txt`, and `cyclotomic_profiles.json`.

Retained canonical fold representatives / folded correlation signatures:

```text
k=25: m=5   17 / 19      m=11   321 / 1020
k=26: m=5   27 / 44      m=11   560 / 1765
k=27: m=5   26 / 45      m=11   660 / 2055
```

## 5. Correlation multiset coverage

For an active weight, the 27 independent correlations are integers with fixed sum and a weight-specific upper bound on their square sum.

**Completeness argument.** The verifier reconstructs all nondecreasing 27-tuples recursively. At a prefix with last value `ell`, remaining length `r`, remaining sum `R`, and accumulated square sum `Q`, it prunes only when

```text
R < ell*r
```

or

```text
Q + balanced_squares(R,r) > C_max.
```

Otherwise it tries every next value from `ell` through `min(k, floor(R/r))`. The manuscript proves by induction on the prefix length that every feasible nondecreasing tuple survives along its own prefix path and is emitted exactly once.

Each emitted multiset represents exactly

```text
27! / product_v multiplicity(v)!
```

ordered profiles.

Before those permutations are materialized, `src/maxdet/boxed_moment.py` supplies a rigorous product bound using an outward dyadic cap on all `q_r`. If the cap proves the entire multiset below the incumbent, all its permutations are discarded at once.

This explains the two profile totals:

```text
coarse ordered count      38,629,684
materialized assignments   6,845,230
unit-canonical profiles       177,988
retained lift targets            2,316
```

The difference between the first two rows is rigorously pruned multiset mass, not an unexplored part of the search.

## 6. Exact formal correlation products

For a symmetric integer correlation profile, the formal paired Fourier factors are algebraic integers. Their product over the 27 conjugate pairs is Galois invariant and therefore an integer even before binary realizability has been established.

`native/order55_correlation_screen.cpp` evaluates that signed integer in two prime fields:

```text
P_1 = 2305843009213696591
P_2 = 2305843009213697141
```

and reconstructs it by signed CRT. The verifier checks, for every materialized multiset, the integer capacity inequality

```text
(2*(55-k))^2 * S2^27 < (P_1*P_2)^2 * 27^27,
```

which implies that the absolute product is strictly below `P_1*P_2/2`. Thus the signed reconstruction is unique.

**Artifacts.** `correlation_profiles.json`, `corr_k*_part*.txt`, and the corresponding target files.

## 7. From 2,316 targets to 16,084 tasks

Each target fixes a full correlation profile. Its residue-class sums select all retained 5-fold and 11-fold vectors with those folded correlations. The verifier explicitly reconstructs those fold-vector sets and takes every compatible pair. Independent row and column translations are normalized through the Chinese-remainder identification

```text
Z/55Z  ~=  Z/5Z x Z/11Z.
```

The resulting target list expands to exactly 16,084 distinct task descriptors. `lift_targets.json` records the expansion and `lift_part*.txt` stores the tasks.

## 8. Complete binary fibers

`native/order55_torus_lift.cpp` writes the word as eleven five-bit column masks. A task fixes five row margins, eleven column margins and the complete periodic autocorrelation.

The eleven columns are split 5+6. A half-word is keyed by the seven additive quantities

```text
five partial row counts + correlations at shifts 11 and 22.
```

A prefix is pruned only when one of those nonnegative partial sums already exceeds its final target. Consequently no prefix of a valid word can be pruned. The right half requests the exact complementary signature. Every joined word is then checked at all remaining independent shifts by 55-bit rotation and population count.

This is the mathematical completeness argument for the binary lift. The optional 6+5 replay is only an implementation consistency check.

**Artifacts.** `lift_targets.json`, `lift_coverage.json`, `lift_part*_audit.jsonl`, and `lift_part*_words.txt`.

Recorded total:

```text
joined words = 15,487,882,832
```

Two normalized witnesses survive, in tasks `2315_0_1` and `2315_1_0`; both complement to weight-28 maximizers and belong to the same affine class.

## 9. Final exact checks

The winning determinant is checked in three representations:

1. finite-field Fourier determinant in `P_1` and `P_2`, followed by signed CRT;
2. an exact integer determinant of the 55×55 circulant matrix using SymPy's domain-aware elimination;
3. the exact cyclotomic factorization `28 * N_5 * N_11 * N_55`.

`winner.json` also stores the complete correlation vector, fold vectors, orbit size and stabilizer size. The affine group has order `55*phi(55)=2200`; direct enumeration gives 2,200 distinct images of the canonical word, so the stabilizer is trivial.

## 10. Certificate integrity

`global.json` has status `COMPLETE` and binds the proof data to `hash_manifest.json`. The recorded SHA-256 of the manifest is

```text
b73a12c5b8502ac5cba20f5164ff351ba5df1a208d7cd2e67f45c359d7eedfc5
```

`source_hashes.json` additionally binds the stored proof data to the exact proof-critical source files. For this reason the journal-oriented documentation revision intentionally does **not** rewrite the already-certified generator/verifier/native source bytes or rename existing `*_audit.json` artifacts. Those historical names mean “run summary/check record”; they do not encode an unresolved logical condition. Rewriting them without regenerating the full certificate would destroy the existing source-hash chain.

`build.json` contains absolute paths from the machine on which the certificate was generated. These paths are provenance only. The verifier identifies the committed native source by basename and checks its SHA-256 before recompilation; the paths themselves do not enter any mathematical comparison.

## 11. Minimal verification commands

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python scripts/verify_order55_global.py
python -m pytest
```

Regenerate the finite enumeration with:

```bash
python scripts/build_order55_global.py --workers 8
python scripts/verify_order55_global.py
```

The manuscript should be read together with this proof map when inspecting the computer-assisted part of the argument.