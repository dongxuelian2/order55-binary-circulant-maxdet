# Exact order-55 circulant maxima v1.0.0

The joint paper by Qichao Wang proves both exact order-55 circulant maxima:

- Binary: `D_01(55) = 134694094094758395331307111329132` (OEIS A086432).
- Sign: `D_pm(55) = 297532404423965431213849494795385191907933028352` (OEIS A215723).
- Normalized sign: `D_pm(55)/2^54 = 16516366298178510386365752382353` (OEIS A215897).

This release contains the single-author joint manuscript and the original
binary and sign proof packages. Verify the stored certificates with
`python scripts/verify_order55_global.py` and
`python scripts/verify_order55_sign.py`. The sign package is stored in Git;
no recompression or certificate regeneration is required. The attached
`SHA256SUMS.txt` gives the PDF digest.
