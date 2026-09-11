# Public-source check, 2026-09-11

- https://oeis.org/A215723 : maximum determinant of circulant (1,-1) matrices;
  linked b-file described as n=1..52.
- https://oeis.org/A215897 : A215723(n)/2^(n-1); linked b-file described as
  n=1..52. Direct b-file fetches returned tool errors, so contents unverified.
- https://arxiv.org/abs/1801.00399 : latest v6, revised 2021-02-20.
- https://arxiv.org/html/1801.00399v6 : Table 4 gives sign normalized maxima
  THROUGH 53, including 512364770126478307153560491081 at n=53.
  This updates the prompt's expected literature frontier.
- https://cs.uwaterloo.ca/journals/JIS/VOL21/Brent/brent11.html :
  original JIS article 18.5.6 (2018).
- General web queries: A215723 55; A215897 55; order 55 circulant maximal
  determinant; maximum determinant circulant sign matrix 55; 55x55 circulant
  determinant. No prior exact sign-55 value identified in returned results.
- Targeted GitHub and arXiv searches likewise found no sign-55 theorem.
  https://github.com/ruu413/aletheia-54 concerns binary order 54 and is outside
  this task. General results about D-optimal order 110 use paired 55-circulants
  and do not solve the single-circulant objective here.

These are bounded search observations, not an absolute-priority claim.
Exact integers and winner words must be searched again before any RESULT.

Repository baseline: 5d3b8e4. Branch sign-circulant-order55, separate worktree
E:/maxdet-sign55. No AGENTS.md found by rg in the original repository.
The original independent_lift_report.json explicitly records
global_completeness=INCOMPLETE and independent_full_fiber_union=NOT_COMPLETED.

