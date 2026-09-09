# Software and provenance audit

## Frozen revision and file boundaries

The required production revision is
`e77a1aa6afea16fe59a35f9c0de93f385f9b8781`.  The audit did not modify
production source, certificate, result, or manuscript files.  All new code and
machine-readable outputs are under `audit/order55/`.  Unrelated tracked edits
and untracked paths already present in the checkout were not staged or edited.

The production certificate reports SHA-256 hashes for 71 certificate
components and seven production source files.  The independent Python checker
recomputed all of those hashes and independently checked the two CRT primes,
the global maximum, and the component cross-links.

## Independent implementations

* `verify_winner_independent.py` is a no-import clean-room Python checker.  It
  constructs the integer circulant directly, runs fraction-free Bareiss, runs
  modular Gaussian elimination at three small primes, evaluates a Fourier
  product at those primes, and computes three cyclotomic resultants.
* `verify_global_independent.py` rederives bounds, fold maps, profile universes,
  CRT products, task topology, and one million random forward maps without
  importing production modules.
* `independent_profile_screen.cpp` is a separate exhaustive native screen with
  a different permutation traversal and independent fold/canonical/CRT code.
* `independent_torus_lift.cpp` uses a typed seven-entry signature and an
  unordered map, rather than the production packed integer/sorted-table path.
  It checks all 27 independent shifts after every join.
* `small_order_ground_truth.cpp` uses Gray-code enumeration and two unrelated
  one-billion-scale primes for exact centered CRT ground truth at orders 15 and
  21.

## Static and dynamic checks

The audit-only C++ files compile with no diagnostics under
`-Wall -Wextra -Wconversion -Wsign-conversion`.  The production native files
were also compiled in syntax-only warning mode.  Their warnings are implicit
signedness or size conversions; inspection of the relevant ranges found no
overflow or undefined behavior, but these should be made explicit.

AddressSanitizer and UndefinedBehaviorSanitizer were available on Windows.
They were run on the audit ground-truth solver at order 15 and on the
clean-room lift solver for a committed order-55 task; both exited successfully
with no sanitizer report.

The repository test suite reports 30 passed tests.  `pip check` reports no
broken requirements.

## Provenance limitations

The certificate does not record a compiler version, native executable hash,
Python executable/environment hash, or a complete project-configuration hash.
The manifest omits `global.json` itself, although `global.json` cross-links the
manifest digest.  Some historical incumbent records carry an earlier
incumbent-search commit; those records are not the final global proof.

These are reproducibility metadata gaps, not evidence of a wrong result.  The
clean-room exhaustive profile and lift replays substantially reduce the risk,
but a publication-quality release should record the exact toolchain, command,
source commit, and executable digest.

## Severity assessment

No unsafe arithmetic was found in the audited ranges and no audit-only
implementation warning remains.  The production conversion warnings and
missing toolchain provenance are MINOR issues.  They do not change the
mathematical verdict.
