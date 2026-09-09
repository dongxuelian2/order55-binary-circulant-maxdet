AUDIT PASS WITH MINOR ISSUES

# Adversarial independent audit: order 55

## Scope and frozen input

The audited claim is

\[
D_{01}(55)=134694094094758395331307111329132,
\]

with one affine-equivalence class of maximizers and 2200 maximizing words.  The
claimed canonical high-weight word is

`0000000100011111011011001101011011010100011110101000111`

of weight 28.  Its complement has weight 27 and determinant
`129883590734231309783760428781663`.

The frozen production revision was
`e77a1aa6afea16fe59a35f9c0de93f385f9b8781`.  The production tracked files
were not edited by this audit.  The current checkout also contains unrelated
tracked edits and untracked work from another task; all of those paths were
left untouched.

## Verdict

The numerical theorem and the uniqueness statement survive an adversarial
recheck.  No fatal mathematical error, omitted active weight, unexplained
profile partition, missing lift task, or false determinant equality was found.
The result is therefore not invalid and does not require a mathematical
revision.

The minor issues are disclosure and provenance issues:

1. The production verifier's advertised “independent” lift check reuses the
   production native lift implementation; changing the split from 5+6 to 6+5
   is a parameter variation, not an independent implementation.
2. The certificate hashes its certificate components and seven source files,
   but does not bind the compiler version, native executable hash, Python
   environment, or all project configuration.  `global.json` is cross-linked
   to the manifest but is not itself a manifest entry.
3. Strict warnings on the production native sources are implicit signedness or
   size conversions.  Range analysis found no resulting overflow or undefined
   behavior, but the conversions should be made explicit before publication.
   There is also a harmless trailing blank line in `src/maxdet/boxed_moment.py`.

These issues affect reproducibility claims and code hygiene, not the audited
maximum.  The first issue is neutralized for this audit by the clean-room
profile and lift implementations listed below.

## Independent evidence

* `verify_winner_independent.py` uses no `maxdet` import.  It agrees on the
  winner by fraction-free Bareiss elimination, modular Gaussian elimination at
  three unrelated small primes, a Fourier product at those primes, and three
  cyclotomic resultant norms.  It checks the complement identity, 2200-word
  affine orbit, stabilizer 1, and 100 deterministic random words.
* `independent_profile_screen.cpp` enumerates every distinct permutation of
  every exhaustive sorted profile using `std::next_permutation`.  It performs
  its own fold tests, unit-canonical test, two-prime signed CRT product, and
  target serialization.  All eight partitions reproduce their committed
  traversal counts and target sets; the total is 2316 targets.
* `independent_torus_lift.cpp` uses a typed seven-coordinate
  `unordered_map` key rather than the production packed/sorted table.  It
  checks all shifts 1 through 27 after the join, including 11 and 22.  A full
  split-6 replay reproduces all 16,084 tasks, 15,487,882,832 joined words, and
  exactly two solutions.
* `small_order_ground_truth.cpp` and `small_order_regression.py` exhaust all
  `2^15=32768` and `2^21=2097152` words, respectively, and independently
  check determinants, complement identities, affine closure, folded profiles,
  and torus fibres.  Both regressions pass.
* The original production verifier was also run with `--full-lifts` and
  exited successfully.  It is recorded as a production self-check, not as
  the independent evidence above.

## Final classification

* Fatal: none found.
* Major revision: none required for the theorem or exhaustive result.
* Minor revision: correct the independence wording, strengthen executable and
  environment provenance, and make production conversion warnings explicit.

The detailed derivation, coverage map, count reconstruction, software review,
and manuscript review are in the neighboring files in this directory.
