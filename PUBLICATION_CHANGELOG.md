# Publication change log

Source sign artifacts: commits `69ec2ed8c4c9b572a37d9d2419247341883891ef`
and `db02b4ccfe409a89feafa2bd984bfd9740d9a7f8` on the preserved
`sign-circulant-order55` branch. The publication branch starts at
`7377ad14e88172822397bb54d3c2c1040780a3cd` and retains those original
objects. No certificate file was recreated or restored from prose.

Publication-facing sources changed: `README.md`, `CITATION.cff`, `paper/main.tex`,
`PROOF_MAP.md`, `docs/OEIS_SUBMISSION.md`, `RELEASE_NOTES.md`, `SHA256SUMS.txt`,
`paper/order55_global/arxiv/main.tex`, `paper/order55_global/manuscript.tex`,
`paper/order55_global/manuscript.md`, all four existing OEIS draft files,
`scripts/check_order55_manuscript.py`, `scripts/check_order55_sign_manuscript.py`,
and `scripts/render_order55_paper.py`. Generated publication PDFs and the
versioned `dist` PDF were rebuilt from the same LaTeX source.

During integration of the newer public main, its deletions of the historical binary audit package were rejected. Original bytes from the sign branch were retained for `certificates/order55_global_v2/`, `certificates/order55_fixed_weight_1_7.*`, `audit/order55/`, `audit/completeness/`, the older binary native and Python utilities, and the binary test files. No proof-critical certificate bytes were edited.
The inherited `.gitattributes` is unchanged.

The old tracked arXiv review bundle was removed by the newer main history; its
Git objects remain recoverable. The historical research script
esearch/order55_sign/upgrade_paper.py retains its previous draft string as
provenance and is not used as current publication metadata.

The v1.0.0 release PDF digest was verified against the commit-pinned
output/pdf/order55_joint.pdf object. OEIS and citation links use that
commit-pinned URL; the release page remains linked for discovery.
