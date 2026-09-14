# Literature audit

## Directly verified primary sources

- Guillermo Nuñez Ponasso, *Combinatorics of Complex Maximal Determinant
  Matrices*, PhD dissertation, Worcester Polytechnic Institute, 2023;
  arXiv:2404.09040v2 (24 February 2026). Section 5.3.3 reports Adam Zsolt
  Wagner's greedy search and prints the order-15 exponent matrix `M15` on
  thesis page 152 / physical PDF page 167. It
  poses its maximality as Research Problem 20. The matrix is transcribed in
  `data/literature_matrices/nunez_ponasso_m15.json`.
- Guillermo Nuñez Ponasso, *Maximal determinants of matrices over the roots of
  unity*, Linear Algebra and its Applications 723 (2025), 201--243,
  DOI 10.1016/j.laa.2025.05.024; arXiv:2503.11114. Table 1 lists, for order
  15, the unproved record `|det|^2/3^14 = 2^22 3^6 19`; the question marks in
  the table explicitly mean maximality is not proved. The paper states that
  exact maxima are calculated only for orders below 14.

## Current interpretation

Exact local verification of the printed matrix gives
`det(M15)=604661760+241864704 omega = 2^11 3^10(5+2 omega)` and hence the
claimed norm. The strongest directly verified published lower bound is therefore
`|det M15|^2 = 2^22 3^20 19`, equivalently
`|det M15| = 2^11 3^10 sqrt(19)`. No later improvement or proof of maximality
has yet been located through 12 September 2026. The February 2026 thesis
revision still states this as Research Problem 20. This is a bounded
literature-search status, not an absolute priority claim. No associated public
code or supplementary data repository was located; the printed matrix is the
accessible primary source.

## Other relevant source

- Arne Winterhof, *On the non-existence of generalized Hadamard matrices*,
  Journal of Statistical Planning and Inference 84 (2000), 337--342,
  DOI 10.1016/S0378-3758(99)00147-0. The 2025 paper cites this for the
  nonexistence of `BH(15,3)`; for this order the short inert-prime norm proof
  in `mathematical_notes.md` is self-contained.
