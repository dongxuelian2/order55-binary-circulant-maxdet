from pathlib import Path
import json,re
ROOT=Path(".");OUT=ROOT/"certificates/order55_sign";PAPER=ROOT/"paper/order55_global"
assert json.loads((OUT/"global.json").read_text())["status"]=="COMPLETE"
assert json.loads((OUT/"binary_companion/full_fiber_report.json").read_text())["status"]=="PASS"
w=json.loads((OUT/"winner.json").read_text());parts=json.loads((OUT/"profile_manifest.json").read_text())["partitions"]
targets=json.loads((OUT/"lift_targets.json").read_text())
counts=[]
for k in range(22,27):
 p=[r for r in parts if r["k"]==k];t=[r for r in targets if r["target"]["k"]==k]
 counts.append((k,sum(r["ordered_count"] for r in p),sum(r["ordered_count"] for r in p if not r["excluded"]),len(t),sum(r["tasks"] for r in t)))
table="\n".join(" & ".join(f"{x:,}" for x in row)+r"\\" for row in counts)
section=r"""
\section{The companion sign-circulant maximum}\label{sec:sign}

Write $D_{\pm}^{\mathrm{circ}}(55)$ for the maximum absolute determinant
over sign circulants and $S_{\pm}(55)=D_{\pm}^{\mathrm{circ}}(55)/2^{54}$.
These are OEIS A215723 and A215897 \cite{oeis-sign,oeis-sign-scaled}.
All sign screening is regenerated; the old binary target set is not
a completeness certificate for this different objective.

\subsection{Exact bridge and complement reduction}
\begin{lemma}[Binary-to-sign bridge]\label{lem:sign-bridge}
For $a\in\{0,1\}^{55}$ of weight $0<k<55$ and $x=2a-\mathbf1$,
\[
 \det C(x)=2^{54}\frac{2k-55}{k}\det C(a),\qquad
 S_{\pm}(55)=\max_{1\le k\le27}\frac{55-2k}{k}D_{\mathrm{circ}}(55,k).
\]
\end{lemma}
\begin{proof}
The identity $C(x)=2C(a)-J$ changes the eigenvalue on the all-ones
line to $2k-55$ and doubles every nontrivial Fourier eigenvalue.
It also holds when the circulant is singular. Conjugate pairing gives
$\det C(a)=k\prod q_r$. Complementing $a$ negates $x$ and preserves
absolute determinant. The constant endpoints have rank one.
\end{proof}
Thus the reduced objective is
\[
 R(a)=(55-2k)\prod_{r=1}^{27}q_r=(55-2k)N_5N_{11}N_{55}.
\]
For the encoding $x=2a-\mathbf1$, the signed row sum is $2k-55$;
$55-2k$ is its magnitude on the reduced weights.
In contrast, the complementary binary objective has coefficient $55-k$.
The small-weight certificates in Appendix A transfer exactly by the
bridge. The binary global theorem and complement duality further give
the exact weight-$27$ sign value
$M/28=4810503360527085547546682547469$.
The sign upper-bound proof below does not require the binary global theorem.

\subsection{Sign correlations and exact bounds}
Put $\rho_s=\sum_i x_ix_{i+s}$. Expansion and orthogonality give
\[
 \rho_s=55-4k+4c_s,\quad \rho_s\equiv3\pmod4,\quad
 55+2\sum_{s=1}^{27}\rho_s=(2k-55)^2.
\]
For $p_r=|f_x(\zeta^r)|^2=4q_r$ and $\sigma=2k-55$,
\[
 \sigma^2+2\sum p_r=55^2,\qquad
 \sigma^4+2\sum p_r^2=55\left(55^2+2\sum_{s=1}^{27}\rho_s^2\right).
\]
At fixed weight this is an invertible affine change of the integer
correlations; the sign congruence adds no condition beyond that
integrality. The implementation keeps the sparse binary deviations.

The first-moment bound is
\[
 R(a)\le(55-2k)\left(\frac{k(55-k)}{54}\right)^{27}.
\]
The fourth-moment and capped-product lemmas above apply to $\prod q_r$,
with every threshold multiplied by $55-2k$. All exclusions are strict
integer or rational outward comparisons against the exact incumbent.
For a mod-$5$/mod-$11$ fold pair $(r,s)$, put
\[
 T_{55}=\frac{55k+k^2-5\sum r_i^2-11\sum s_i^2}{2}.
\]
The joint-fold bound is
\[
 (55-2k)N_5(r)N_{11}(s)(T_{55}/20)^{20}.
\]
The primitive fourth moment subtracts the two folded fourth moments
from the full fourth moment, so the same product lemma can be used
in that sector. The recorded joint screen combines the first-sector
bound with the global fourth-moment bound; all conservative survivors
are subsequently resolved by complete profiles and lifts.

\subsection{Complete screening and independent audit}
The exact moment screen leaves $k=22,23,24,25,26$. Capped moments
reduce $9,531,797,706$ ordered formal correlations to $731,378,844$
explicit evaluations in $17$ multisets.
\begin{center}
\small
\begin{tabular}{@{}rrrrr@{}}
\toprule
$k$ & Formal profiles & Explicit profiles & Targets & Lift tasks\\
\midrule
TABLE
\bottomrule
\end{tabular}
\end{center}
There are $17,835$ canonical targets and $151,772$ lift tasks.
No target remains at weight $22$.

The independent arithmetic checker parametrizes the two stationary
values by their separation, verifies its own dyadic cosine enclosures,
and generates multisets by descending multiplicities.
An independent bounded-alphabet FKM necklace generator produces every
cyclic fold representative; it shares no production affine-minimum test.
All fold representative sets and their correlations agree exactly.

An independent \texttt{next\_permutation} traversal computes every
canonical algebraic product by CRT without production interval-product
pruning. Its raw above-threshold set contains $255$ additional indefinite
formal spectra. For every one, a recorded frequency and exact rational
upper bound prove $q_j<0$, excluding a binary autocorrelation.
After these explicit positive-semidefinite discharges, the independent
and production target sets, including their products, are exactly equal.

The production five-plus-six lift uses a sorted packed integer key.
The independent six-plus-five implementation uses typed seven-coordinate
hash keys with equality checks and tests every shift. Both exhaust all
$151,772$ tasks and inspect $113,997,335,580$ joined words.
Their joined counts and complete output sets agree task by task,
returning exactly two normalized words.

For global coverage, map any feasible word's correlation to its
canonical unit representative. Its two folds then occur in the
independent necklace lists in that orientation. The two fold rotations
are simultaneously induced by one translation modulo $55$, by CRT.
Consequently the Cartesian products of listed margin fibers contain
a transformed representative of every feasible word. An independent
task audit reconstructs all $151,772$ descriptors from those fold lists
and checks their exact partition. The two full lift enumerations
therefore prove global exhaustiveness.

\subsection{Winner and all equality cases}
Both returned words are affine-equivalent to the weight-$23$ word
\begin{lstlisting}
0000000000100110110100111001001110001010101100010010111
\end{lstlisting}
whose sign encoding is
\begin{lstlisting}
----------+--++-++-+--+++--+--+++---+-+-+-++---+--+-+++
\end{lstlisting}
Its exact norms are
\[
 N_5=131,\quad N_{11}=917467,\quad N_{55}=15268987821561877723321,
\]
and its binary determinant is
\[
 \det C(a)=42208491650900637654045811643791.
\]
Direct integer Bareiss on the sign matrix gives a negative determinant
with absolute value
\[
 297532404423965431213849494795385191907933028352.
\]
Modular Gaussian elimination over three independently checked primes
and signed CRT agree, with sufficient capacity proved by Hadamard's
bound. Both agree with $2^{54}(2k-55)N_5N_{11}N_{55}$.

An independent pushforward-support orbit calculation gives $2200$
words and trivial stabilizer. Global negatives have weight $32$ and
row sum $9$, forming a disjoint second affine orbit, with canonical
binary representative
\begin{lstlisting}
0000001101110110101010100110111101111011110011010001111
\end{lstlisting}
There are exactly $4400$ maximizing sign words, two affine orbits
of size $2200$ with stabilizer one, and one affine-plus-negation class.
This proves Theorem~\ref{thm:sign} and classifies all equality cases.

\subsection{Comparison and regression}
The sign winner is not induced by the binary global winner.
Its weight-$23$ binary determinant falls below $D_{01}(55)$ by
$92485602443857757677261299685341$. Its weight-$32$ complement has
binary determinant $58724857949079148040411564026144$.
The half binary correlations are $9$ seventeen times and $10$ ten times;
the $54$ nontrivial sign correlations are $-1$ thirty-four times and
$3$ twenty times. Relative to the old binary winner, the product of
nontrivial Fourier pairs decreases but the absolute row-sum factor
increases from $1$ to $9$. The sign objective improves by approximately
$3.43$, the ratio of the two exact normalized integers.

At orders $15$ and $21$, direct Gray-code sign Fourier enumeration
covers all $2^{15}$ and $2^{21}$ words. The reduced sign moment/profile
pipeline followed by a CRT-column split lift agrees on the maximum,
weight reduction, profile coverage, lift coverage and complete winner
sets. The normalized maxima are $23859$ and $39337984$, with $480$
and $84$ winning words. A further $707$ exact bridge and moment checks
include all words at $3,5,7,9$ and one sampled word at every reduced
weight at order $55$.

\subsection{Public sources and reproducibility}
The check on 11 September 2026 inspected both OEIS entries, the JIS
article and the latest Brent--Yedidia revision. The OEIS linked tables
are described as ending at $52$, whereas arXiv v6 Table 4 gives sign
values through $53$. Searches of the exact new integers, both encodings,
arXiv, GitHub, Zenodo and the general web located no earlier exact
sign-$55$ determination. To the best of our knowledge this is the
missing companion result; the search does not establish absolute priority.

All sign certificates are under \path|certificates/order55_sign/|.
The ordinary verifier rechecks hashes, exact winner arithmetic, weight
bounds, target sets, negative-spectrum certificates, task records and
orbit classification. The second command recompiles and independently
reruns the complete proof without any heuristic search:
\begin{lstlisting}
python scripts/verify_order55_sign.py
python scripts/verify_order55_sign.py --replay
\end{lstlisting}
"""
section=section.replace("TABLE",table)
texpath=PAPER/"arxiv/main.tex";tex=texpath.read_text()
tex=tex.replace("pdftitle={The Maximum Determinant of a 55 by 55 Binary Circulant Matrix}","pdftitle={Maximum Determinants of 55 by 55 Circulant Matrices over 0,1 and -1,1}")
tex=tex.replace(r"\title{The Maximum Determinant of a $55\times55$ Binary Circulant Matrix}",
r"\title{Maximum Determinants of $55\times55$ Circulant Matrices\\over $\{0,1\}$ and $\{-1,1\}$}")
abstract=r"""\begin{abstract}
We determine the maximum determinants of circulant matrices of order $55$
over both binary alphabets:
\[
 D_{01}(55)=134694094094758395331307111329132,
\]
\[
 D_{\pm}^{\mathrm{circ}}(55)=
 297532404423965431213849494795385191907933028352.
\]
The normalized sign maximum is $16516366298178510386365752382353$.
The binary maximizers have weight $28$ and form one affine orbit of
$2200$ words. The sign maximizers have $23$ or $32$ plus signs and form
two affine orbits of $2200$ words each, paired by global negation.
The exact binary-to-sign bridge shows why these objectives have different
winners. The proofs combine spectral moments, capped product bounds,
cyclotomic folds, correlation-profile enumeration and binary lifts.
The sign screen retains $17,835$ targets and $151,772$ margin tasks.
Production and independent lifts each inspect $113,997,335,580$ joined
words and agree on every task and every returned word. Independent
profile generation, explicit negative-spectrum certificates,
small-order full-domain regressions and three exact determinant
calculations complete the audit. To the best of our knowledge neither
exact order-$55$ maximum was previously published in the sources checked.
\end{abstract}"""
tex=re.sub(r"\\begin\{abstract\}.*?\\end\{abstract\}",lambda m:abstract,tex,flags=re.S)
headline=r"""
\begin{theorem}[Sign maximum and classification]\label{thm:sign}
For sign circulants of order $55$,
\[
 D_{\pm}^{\mathrm{circ}}(55)=
 297532404423965431213849494795385191907933028352,
\]
\[
 S_{\pm}(55)=D_{\pm}^{\mathrm{circ}}(55)/2^{54}
 =16516366298178510386365752382353.
\]
All maximizers have $23$ or $32$ plus signs, and row sum $-9$ or $9$.
They form two affine orbits, each of size $2200$ with trivial stabilizer,
and one class under affine maps and global negation.
\end{theorem}
"""
tex=tex.replace(r"\end{theorem}",r"\end{theorem}"+"\n"+headline,1)
tex=tex.replace("The main result is the following.","We obtain the following two headline results.")
tex=tex.replace("We determine the next order, $55=5\\cdot11$.","We determine both alphabet variants at order $55=5\\cdot11$.")
old="""is a split-dependent cross-check because it reuses the same native lift
algorithm; it is not an independent full fiber-union proof.  A separate small
fiber is cross-checked against all $78,125$ direct column assignments using
three splits.  The clean-room audit therefore certifies the listed task
partition and the displayed witnesses, while leaving the independent union
certificate for every full order-$55$ fiber open."""
new=r"""was initially a split-dependent cross-check. In this joint revision,
the separate typed-key implementation exhausts all $16,084$ original
binary tasks with split six-plus-five. All $15,487,882,832$ joined words
and both normalized solutions agree task by task with production.
The complete supplement is
\path|certificates/order55_sign/binary_companion/full_fiber_report.json|.
A separate small fiber is checked against all $78,125$ direct column
assignments using three splits."""
assert old in tex;tex=tex.replace(old,new)
old="""and have the same affine canonical representative.  The implication from
the finite profile screen to all binary words is conditional on the full
fiber-union statement; the current review certificate does not silently
replace that statement by the split-dependent replay."""
new="""and have the same affine canonical representative. The newly completed
independent full-fiber computation discharges the remaining union
condition in the historical review snapshot."""
assert old in tex;tex=tex.replace(old,new)
start=tex.index(r"\begin{lemma}[Lift completeness, conditional form]");end=tex.index(r"\section{Reproducibility and availability}",start)
tex=tex[:start]+r"""\begin{lemma}[Lift completeness]
Every binary word with determinant at least $M$ has a transformed
representative in the union of the $16,084$ listed margin fibers,
each of which is exhausted by production and independent implementations.
The preceding profile and pruning lemmas therefore imply
Theorem~\ref{thm:main}.
\end{lemma}
\begin{proof}
The target and task sets were independently reconstructed in the
earlier review. The new supplement exhausts their entire fibers with
the separate typed-key implementation and checks exact task counts
and word sets. CRT combines the two fold rotations into a translation
modulo $55$, and unit normalization preserves the determinant.
Thus the listed fibers cover a representative of every eligible word.
\end{proof}

The historical v2 independent lift report is preserved unchanged.
Its previously missing full-fiber step is now supplied by
\path|certificates/order55_sign/binary_companion/full_fiber_report.json|.
The binary value and frozen original certificates are unchanged.

""" +tex[end:]
tex=tex.replace(r"\section{Conclusion}",section+"\n"+r"\section{Conclusion}",1)
tex=tex.replace("with one affine class of $2200$ maximizers.",
"""with one affine class of $2200$ maximizers. The companion sign maximum is
$16516366298178510386365752382353\\cdot2^{54}$, with $4400$ words in
two affine orbits paired by negation.""")
tex=tex.replace("""The source code and machine-readable certificates accompanying this
manuscript will be archived with the public version of this work.  A
permanent repository URL was not available at packaging time, so no
unresolved URL placeholder is included in the arXiv source.""",r"""The project repository is
\url{https://github.com/dongxuelian2/order55-binary-circulant-maxdet}.
The joint revision is prepared on \path|sign-circulant-order55|.
This task does not publish or submit the revision; local OEIS drafts
accompany the source.""")
texpath.write_text(tex)
mdpath=PAPER/"manuscript.md";md=mdpath.read_text()
md=md.replace("# The Maximum Determinant of a 55 × 55 Binary Circulant Matrix","# Maximum Determinants of 55 × 55 Circulant Matrices over {0,1} and {-1,1}",1)
md=md.replace("Reproducible computational research note — 9 September 2026","Qichao Wang and Daoyu Dong; equal contribution. Joint revision, 11 September 2026.",1)
note="""The companion sign maximum is now also exact:
297532404423965431213849494795385191907933028352, or
16516366298178510386365752382353 after division by 2^54.
Its 4,400 maximizing words have weights 23 and 32 and form two affine
orbits paired by negation. All new independent audits pass.

"""
md=md.replace("## 1. Introduction and previous exact frontier",note+"## 1. Introduction and previous exact frontier",1)
md=md.replace("The current package is a review draft, not an arXiv-ready submission: an independent full fiber-union certificate (or equivalent branch-level partition/union proof) is still required. The exact theorem and incumbent are retained without weakening them.",
"The independent full binary fiber-union supplement is now complete in certificates/order55_sign/binary_companion/full_fiber_report.json; the historical v2 review snapshot is preserved unchanged.")
md=md.replace("while the full order-55 fiber-union\ncertificate remains incomplete.","and now includes the completed independent full order-55 fiber-union supplement.")
mdsection=r"""
## 13. Exact sign-circulant companion theorem

**Theorem 2.**
\[
D_{\pm}^{\mathrm{circ}}(55)=297532404423965431213849494795385191907933028352,
\]
\[
S_{\pm}(55)=16516366298178510386365752382353.
\]
The winning weights are 23 and 32, with signed row sums -9 and 9.
There are 4,400 maximizing sign words: two affine orbits of size 2,200
with trivial stabilizers, and one affine-plus-negation class.

For x=2a-1, the trivial Fourier eigenvalue is 2k-55 and all others double:
\[
\det C(2a-\mathbf1)=2^{54}\frac{2k-55}{k}\det C(a).
\]
Consequently
\[
S_{\pm}(55)=\max_{1\le k\le27}\frac{55-2k}{k}D_{\mathrm{circ}}(55,k).
\]
The complete proof uses newly derived sign thresholds throughout.

The moment bounds leave k=22,...,26. Capped bounds reduce 9,531,797,706
formal ordered correlations to 731,378,844 explicit evaluations in
17 multisets. They retain 17,835 canonical targets and 151,772 margin
fibers. Both independent lift implementations inspect 113,997,335,580
joined words and agree on every task and both normalized outputs.

The independent audit reconstructs all multisets, FKM-necklace fold
representatives, and ordered profiles. Its 255 extra indefinite
algebraic spectra are each discharged by an exact negative-eigenvalue
upper bound. The realizability-relevant target sets then match exactly.
An independent reconstruction of every Cartesian margin fiber, together
with simultaneous CRT translation of the two folds, proves coverage of
all binary words and hence all sign words.

Canonical binary and sign encodings:
~~~text
0000000000100110110100111001001110001010101100010010111
----------+--++-++-+--+++--+--+++---+-+-+-++---+--+-+++
~~~
The second affine orbit has canonical binary encoding
~~~text
0000001101110110101010100110111101111011110011010001111
~~~
Bareiss, modular Gaussian elimination plus CRT, and cyclotomic
resultants agree, with
~~~text
N5 = 131
N11 = 917467
N55 = 15268987821561877723321
~~~
The weight-23 binary determinant is 42208491650900637654045811643791,
a shortfall of 92485602443857757677261299685341 from the binary maximum.
Its weight-32 complement has determinant 58724857949079148040411564026144.
Binary nontrivial correlations are 9 thirty-four times and 10 twenty
times; sign correlations are -1 thirty-four times and 3 twenty times.
The absolute row-sum factor grows from 1 to 9, compensating for the
smaller nontrivial Fourier product.

For sigma=2k-55, rho_s=55-4k+4c_s and p_r=4q_r,
\[
55+2\sum\rho_s=\sigma^2,\quad \sigma^2+2\sum p_r=55^2,\quad
\sigma^4+2\sum p_r^2=55(55^2+2\sum\rho_s^2).
\]
All sums here run through the 27 paired nonzero shifts or frequencies.
The sign congruence rho_s=3 mod 4 is equivalent to binary integrality.

Full-domain sign regressions at 15 and 21 agree on maximum, weight
reduction, profile coverage, lift coverage and every winning word:
normalized maxima 23859 and 39337984, with 480 and 84 winners.
The bridge/moment regression passes 707 exact cases.

Verify the certificate or recompile and replay the full proof:
~~~powershell
python scripts/verify_order55_sign.py
python scripts/verify_order55_sign.py --replay
~~~
The formal typeset derivation is in paper/order55_global/arxiv/main.tex.
The public-source check on 11 September 2026 found no prior exact sign-55
value in the inspected sources. OEIS links tables through 52, whereas
Brent-Yedidia arXiv v6 Table 4 includes 53. Exact values and encodings
were searched again. To the best of our knowledge this is the missing
companion result; no absolute-priority claim is made.

"""
md=md.replace("## Appendix A.",mdsection+"## Appendix A.",1);mdpath.write_text(md)
bib=PAPER/"arxiv/references.bib";bib.write_text(bib.read_text()+r"""
@misc{oeis-sign,
 author={{OEIS Foundation Inc.}},
 title={A215723: Maximum determinant of an n by n circulant (1,-1)-matrix},
 year={2026}, note={Inspected 11 September 2026},
 url={https://oeis.org/A215723}
}
@misc{oeis-sign-scaled,
 author={{OEIS Foundation Inc.}},
 title={A215897: A215723(n) divided by 2 to the power n-1},
 year={2026}, note={Inspected 11 September 2026},
 url={https://oeis.org/A215897}
}
""")
(PAPER/"REVIEW_RESPONSE.md").write_text("""# Response to the global-completeness review

Status: RESOLVED in the joint revision of 11 September 2026.

The original binary maximum and frozen binary certificates are preserved.
The typed-key independent lift now exhausts all 16,084 original binary
fibers with split 6+5, agreeing on 15,487,882,832 joined words and both
normalized solutions. The full report is stored under
certificates/order55_sign/binary_companion/full_fiber_report.json.
The old v2 report remains a historical snapshot.

The sign theorem has a wholly regenerated target universe. Its independent
audit checks moment bounds, cosine enclosures, multisets, all FKM-necklace
folds and exact target equality after 255 certified indefinite spectra
are excluded. Its full fiber audit reconstructs all 151,772 tasks and
compares two independent enumerations of 113,997,335,580 joined words.

Full-domain sign regressions at 15 and 21 and 707 bridge/moment checks pass.
Independent pushforward orbit enumeration gives 4,400 maximizing sign
words, two affine orbits and one class with global negation.

No publication or OEIS submission has been performed.
""")
print("JOINT MANUSCRIPT AUTHORED",counts)

