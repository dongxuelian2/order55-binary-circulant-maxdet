"""Exact fourth moment with a rigorous uniform spectral cap from rearrangement."""
from fractions import Fraction as F
from maxdet.order55 import sqrt_interval
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1]
M=103755181499155107923236157710676
cos=[tuple(map(int,line.split())) for line in (ROOT/'results/order55_cosine_intervals.txt').read_text().splitlines()]
S=1<<24

from maxdet.boxed_moment import boxed

out=[]
for k in (24,25,26,27):
    parts=json.loads((ROOT/f'results/order55_correlation_partitions_k{k}.json').read_text())
    for idx,p in enumerate(parts['partitions']):
        base=k*(k-1)//54
        d=sorted(v-base for v in p['values'])
        caps=[]
        for j in (1,5,11):
            # Sort exact cosine values via their disjoint dyadic enclosures.
            cs=sorted(cos[min(j*s%55,55-j*s%55)] for s in range(1,28))
            caps.append(F(k-base)+sum(F(v*(pair[1] if v>=0 else pair[0]),S) for v,pair in zip(d,cs)))
        cap=max(caps)
        s1=F(k*(55-k),2);s2=F(55*(k*k+2*p['squares'])-k**4,2)
        bound=(55-k)*boxed(s1,s2,27,cap)
        row=dict(k=k,partition=idx,cap=str(cap),upper=int(bound),excluded=int(bound)<M)
        out.append(row)
        print(row,flush=True)
(ROOT/'results/order55_boxed_moment_screen.json').write_text(json.dumps(out,indent=2))
