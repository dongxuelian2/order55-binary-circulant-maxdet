# Order-55 arXiv review bundle

Release status: **MAJOR REVISION REQUIRED**.  This is an Outcome B review
bundle, not a submission-ready arXiv package.

The exact theorem claim is retained:

```text
D_01(55) = 134694094094758395331307111329132
```

The claimed maximizers form one affine class of 2,200 weight-28 words, with
canonical representative

```text
0000000100011111011011001101011011010100011110101000111
```

The external completeness review found that the existing five-plus-six lift
and six-plus-five replay share the same native lift implementation.  The
second run is therefore only a split-dependent cross-check, not an independent
full fiber-union proof.  The missing independent order-55 fiber-union
certificate is the sole release blocker recorded here; the theorem has not
been weakened to a candidate statement.

## What is included

- `paper/order55_global.pdf`: the revised 11-page manuscript;
- `arxiv/arxiv_submission.zip`: clean `main.tex` + `references.bib` source;
- `paper/REVIEW_RESPONSE.md`: point-by-point response to the four review concerns;
- `certificates/`: compact final review reports, aggregate gate, and compressed pruning table;
- `audit/`: clean-room source files used by the independent checks;
- `verification/`: verifier source, commands, and captured outputs;
- `provenance/`: production commit, packaging metadata, and SHA-256 manifest.

## Verification status

- production verifier: PASS;
- split-dependent `--full-lifts` replay: PASS;
- clean-room profile target set: PASS;
- listed lift task partition: PASS;
- independent witness and affine audits: PASS;
- n=15 and n=21 brute-force reduction regressions: PASS;
- 1,000,000-word forward coverage sample: PASS, unknown=0;
- independent full order-55 fiber union: NOT COMPLETED.

The production theorem commit is
`e77a1aa6afea16fe59a35f9c0de93f385f9b8781`.  The final packaging commit is reported at handoff and can be recovered with
`git rev-parse HEAD` after checkout.

## Fast checks from the repository root

```powershell
& ./.venv/Scripts/python.exe scripts/verify_order55_global.py
& ./.venv/Scripts/python.exe scripts/verify_structured.py
& ./.venv/Scripts/python.exe -m pytest
& ./.venv/Scripts/python.exe -m pip check
& ./.venv/Scripts/python.exe audit/completeness/independent_lift_checker.py
& ./.venv/Scripts/python.exe audit/completeness/affine_canonicalization_audit_v2.py
& ./.venv/Scripts/python.exe audit/completeness/small_order_validation.py
& ./.venv/Scripts/python.exe audit/completeness/random_forward_coverage_v8.py
```

The full production replay is:

```powershell
& ./.venv/Scripts/python.exe scripts/verify_order55_global.py --full-lifts
```

The clean-room profile evaluator is more expensive and writes a compressed
machine-readable table (`certificates/pruned_profiles_v3.jsonl.gz`); use the
exact command in `verification/commands.txt`.

## Reviewer feedback convention

Please label feedback as `FATAL`, `MAJOR`, `MINOR`, or `TYPO`.  In particular,
the unresolved fiber-union item should be evaluated as a `MAJOR` completeness
issue until a second generator or a branch-level partition/union proof is
supplied.

## arXiv metadata choices still requiring author confirmation

- arXiv account and submission ownership;
- subject classification (suggested primary: `math.CO`);
- license choice;
- whether to disclose the exact computational-assistance wording in the
  comments field;
- final public repository/DOI URL, if one is created.

No email addresses, credentials, private URLs, or internal prompts are included
in this bundle.
