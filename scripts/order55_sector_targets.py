"""Sector-wise fourth moments for surviving correlation targets."""
from fractions import Fraction as F
from functools import lru_cache
from pathlib import Path
import json
from maxdet.order55 import moment_product_upper
ROOT=Path(__file__).resolve().parents[1]
M=103755181499155107923236157710676
moment_product_upper=lru_cache(None)(moment_product_upper)
@lru_cache(None)
def fold_dictionary(k,m):
    lookup={}
    for line in (ROOT/f'results/order55_profiles_m{m}_k{k}.jsonl').read_text(encoding='utf-8-sig').splitlines():
        row=json.loads(line);a=row['profile']
        h=[sum(a[t]*a[(t+s)%m] for t in range(m)) for s in range(m)]
        for u in range(1,m):
            key=tuple(h[u*s%m] for s in range(m))
            if key in lookup: assert lookup[key]==row['norm']
            lookup[key]=row['norm']
    return lookup

def apply(k,path):
    output=[];count=0;best=0
    for line in path.read_text(encoding='utf-8-sig').splitlines():
        if not line:continue
        t=json.loads(line);count+=1
        c=[k]+t['correlations']+list(reversed(t['correlations']))
        s1=F(k*(55-k),2);s2=F(55*sum(x*x for x in c)-k**4,2)
        norm=1
        for m in (5,11):
            h=tuple(sum(c[i::m]) for i in range(m))
            norm*=fold_dictionary(k,m)[h]
            s1-=F(m*h[0]-k*k,2)
            s2-=F(m*sum(x*x for x in h)-k**4,2)
        bound=(55-k)*norm*moment_product_upper(s1,s2,20)
        best=max(best,int(bound))
        if int(bound)>=M:
            t['primitive_S1']=str(s1);t['primitive_S2']=str(s2);t['sector_upper']=str(int(bound));output.append(t)
    dest=path.with_name(path.stem+'_sector.jsonl')
    dest.write_text(''.join(json.dumps(t)+'\n' for t in output))
    print(path.name,'input',count,'survivors',len(output),'upper',best,flush=True)
    return dict(path=str(path),input=count,survivors=len(output),upper=best)

if __name__=='__main__':
    import sys
    k=int(sys.argv[1])
    if len(sys.argv)>2: paths=[Path(sys.argv[2])]
    elif k==25:paths=[ROOT/'results/order55_correlation_targets_k25.jsonl']
    else:paths=sorted((ROOT/'results/order55_correlation_work').glob(f'k{k}_part??_targets.jsonl'))
    for p in paths:apply(k,p)
