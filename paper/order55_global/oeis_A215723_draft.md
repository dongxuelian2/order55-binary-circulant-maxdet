# OEIS A215723: prepared update draft

Status: not submitted. The proof and joint manuscript are complete locally.
Publish the prepared branch/certificate revision before sending these links
to OEIS; this task did not push the repository or submit an update.

## Proposed exact term

~~~text
55 297532404423965431213849494795385191907933028352
~~~

Definition: Maximum absolute determinant of a 55 by 55 circulant {-1,+1} matrix.

## Proposed comment

a(55) = 297532404423965431213849494795385191907933028352. The maximum is attained by the circulant whose first
row has binary encoding
0000000000100110110100111001001110001010101100010010111
with 0 interpreted as -1 and 1 as +1. Every maximizing sign word has
23 or 32 plus signs. There are exactly 4400 maximizing first rows, in
two affine orbits of size 2200 with trivial stabilizers, interchanged
by global sign reversal.

The result is a global exact maximum over all sign circulants of order 55.
The proof regenerates all objective-dependent correlation bounds, excludes
every remaining nonwinning profile by full binary lifts, and includes
independent profile and full fiber-union audits.

## Attribution and references

Qichao Wang and Daoyu Dong, equal contribution, 11 September 2026.
Supporting paper: Maximum Determinants of 55 x 55 Circulant Matrices
over {0,1} and {-1,1}.

- Joint paper source (prepared branch):
  https://github.com/dongxuelian2/order55-binary-circulant-maxdet/blob/sign-circulant-order55/paper/order55_global/arxiv/main.tex
- Certificate revision: db02b4ccfe409a89feafa2bd984bfd9740d9a7f8
- Exact certificate:
  https://github.com/dongxuelian2/order55-binary-circulant-maxdet/blob/db02b4ccfe409a89feafa2bd984bfd9740d9a7f8/certificates/order55_sign/global.json
- Verifier:
  https://github.com/dongxuelian2/order55-binary-circulant-maxdet/blob/db02b4ccfe409a89feafa2bd984bfd9740d9a7f8/scripts/verify_order55_sign.py
- Brent and Yedidia, JIS 21 (2018), Article 18.5.6, updated arXiv v6:
  https://arxiv.org/abs/1801.00399v6

## Exact normalization and scope

~~~text
A215897(55) = 16516366298178510386365752382353
A215723(55) = 297532404423965431213849494795385191907933028352
A215723(55) = 2^54 * A215897(55)
~~~

This draft proposes only the n=55 term and its evidence. Existing terms
and the separate A086432 draft are preserved. No missing intermediate
term is inferred from this result.
