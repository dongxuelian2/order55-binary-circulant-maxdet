# Order-55 proof map

| Result | Stored certificate | Normal verifier | Optional full replay |
| --- | --- | --- | --- |
| Binary `134694094094758395331307111329132` | `certificates/order55_global/` | `python scripts/verify_order55_global.py` | `python scripts/verify_order55_global.py --full-lifts` |
| Sign `297532404423965431213849494795385191907933028352` | `certificates/order55_sign/` | `python scripts/verify_order55_sign.py` | `python scripts/verify_order55_sign.py --replay` |

The sign normalized value is `16516366298178510386365752382353`. The sign certificate has a
SHA-256 asset manifest and source-hash manifest, both checked by its normal
verifier. Its production and independent lift reports both record 151,772
tasks and 113,997,335,580 joined words. The binary companion full-fiber
report is `certificates/order55_sign/binary_companion/full_fiber_report.json`.

Original sign certificate and verifier content was committed at
`69ec2ed8c4c9b572a37d9d2419247341883891ef`; portable hash preservation
was committed at `db02b4ccfe409a89feafa2bd984bfd9740d9a7f8`.
The binary production theorem was committed at
`e77a1aa6afea16fe59a35f9c0de93f385f9b8781`.
