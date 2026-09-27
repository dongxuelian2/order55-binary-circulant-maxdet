# Exact order-55 circulant determinants

The exact maxima for order-55 circulant matrices are:

| Alphabet | Exact maximum | OEIS |
| --- | ---: | --- |
| Binary {0,1} | `D_01(55) = 134694094094758395331307111329132` | [A086432](https://oeis.org/A086432) |
| Sign {-1,1} | `D_pm(55) = 297532404423965431213849494795385191907933028352` | [A215723](https://oeis.org/A215723) |
| Sign, divided by `2^54` | `16516366298178510386365752382353` | [A215897](https://oeis.org/A215897) |

The binary maximizers have weight 28: one affine orbit of 2,200 words.
The sign maximizers have weights 23 and 32: 4,400 words in two affine
orbits paired by global negation.

## Publication and proof packages

- Joint single-author manuscript: [LaTeX source](paper/order55_global/arxiv/main.tex)
  and [commit-pinned PDF](https://raw.githubusercontent.com/dongxuelian2/order55-binary-circulant-maxdet/e2204cc7d784cb264f47f0820b1fb71487cb0ecc/output/pdf/order55_joint.pdf) and [v1.0.0 release](https://github.com/dongxuelian2/order55-binary-circulant-maxdet/releases/tag/order55-v1.0.0).
- Binary proof: [certificate](certificates/order55_global/global.json),
  [verifier](scripts/verify_order55_global.py), and [proof map](PROOF_MAP.md).
- Sign proof: [certificate](certificates/order55_sign/global.json),
  [hash manifest](certificates/order55_sign/hash_manifest.json),
  [verifier](scripts/verify_order55_sign.py), and
  [independent full-fiber report](certificates/order55_sign/cleanroom_lift_report.json).
- The [binary full-fiber supplement](certificates/order55_sign/binary_companion/full_fiber_report.json)
  completes the historical v2 review snapshot. The original binary
  certificates remain byte-for-byte unchanged.
- [OEIS submission text](docs/OEIS_SUBMISSION.md) gives comments and links
  for all three entries. It proposes no missing intermediate terms.

The sign certificate records 17,835 correlation targets and 151,772 lift
tasks. The production and independent lifts each examined 113,997,335,580
joined words and agreed task by task. These are stored replay reports; the
normal verifier checks their hashes, structure, and exact arithmetic.

## Verify

Use Python 3.11+ with the dependencies in `pyproject.toml`:

```text
python scripts/verify_order55_global.py
python scripts/verify_order55_sign.py
python -m pytest
python -m pip check
```

The optional full replays are computationally expensive:

```text
python scripts/verify_order55_global.py --full-lifts
python scripts/verify_order55_sign.py --replay
```

The sign package was preserved from Git commits `69ec2ed8c4c9b572a37d9d2419247341883891ef`
and `db02b4ccfe409a89feafa2bd984bfd9740d9a7f8`. The binary theorem
originates at `e77a1aa6afea16fe59a35f9c0de93f385f9b8781`.
