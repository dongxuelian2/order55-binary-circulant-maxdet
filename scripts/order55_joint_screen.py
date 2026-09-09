"""Combine both exact folded norms and integer correlation orbit sums."""
from fractions import Fraction
from math import factorial
from functools import lru_cache
import json
from pathlib import Path
from maxdet.order55 import balanced_squares, fourth_bound
ROOT=Path(__file__).resolve().parents[1]
screen=json.loads((ROOT/'results/order55_weight_screen.json').read_text())
M=int(screen['incumbent'])
fourth_bound=lru_cache(None)(fourth_bound)
out=[]
for k in range(24,28):
    p5=[json.loads(x) for x in (ROOT/f'results/order55_profiles_m5_k{k}.jsonl').read_text(encoding='utf-8-sig').splitlines()]
    p11=[json.loads(x) for x in (ROOT/f'results/order55_profiles_m11_k{k}.jsonl').read_text(encoding='utf-8-sig').splitlines()]
    survivors=[]
    best=0
    for r in p5:
        for s in p11:
            A=Fraction(r['squares']-k,2) # five +/- paired shifts divisible by 5
            B=Fraction(s['squares']-k,2) # two pairs divisible by 11
            C=Fraction(k*(k-1),2)-A-B
            if min(A,B,C)<0 or any(x.denominator!=1 for x in (A,B,C)):
                continue
            minhalf=balanced_squares(int(A),5)+balanced_squares(int(B),2)+balanced_squares(int(C),20)
            # First moments in order-5 and order-11 frequency sectors.
            T=Fraction(55*k+k*k-5*r['squares']-11*s['squares'],2)
            if T<=0:continue
            bound=(55-k)*r['norm']*s['norm']*(T/20)**20
            bound=min(bound,fourth_bound(k,minhalf))
            best=max(best,int(bound))
            if int(bound)<M:continue
            survivors.append(dict(mod5=r['profile'],mod11=s['profile'],N5=r['norm'],N11=s['norm'],
                                  primitive_first_moment=str(T),orbit_correlation_sums=[int(A),int(B),int(C)],
                                  minimum_half_squares=minhalf,upper_floor=int(bound)))
    (ROOT/f'results/order55_joint_k{k}.json').write_text(json.dumps(survivors))
    row=dict(k=k,profile_pairs=len(p5)*len(p11),survivors=len(survivors),upper_floor=best)
    out.append(row)
    print(row,flush=True)
(ROOT/'results/order55_joint_summary.json').write_text(json.dumps(dict(incumbent=str(M),status='INCOMPLETE: no binary lift coverage',weights=out),indent=2))

# Enumerate multisets of paired correlations; count ordered profiles exactly.
for k in range(24,28):
    entry=screen['weights'][k-1]
    limit=entry['maximum_half_correlation_squares_not_excluded']
    total=k*(k-1)//2
    partitions=[]
    def walk(lo,left,remaining,squares,values):
        if not left:
            if remaining==0:
                from collections import Counter
                count=factorial(27)
                for c in Counter(values).values():count//=factorial(c)
                partitions.append(dict(values=values,squares=squares,ordered_count=count))
            return
        if remaining<lo*left or squares+balanced_squares(remaining,left)>limit:return
        for v in range(lo,min(k,remaining//left)+1):
            walk(v,left-1,remaining-v,squares+v*v,values+[v])
    walk(0,27,total,0,[])
    payload=dict(k=k,maximum_half_squares=limit,partitions=partitions,ordered_count=sum(p['ordered_count'] for p in partitions),status='multiset enumeration only; realizability unresolved')
    (ROOT/f'results/order55_correlation_partitions_k{k}.json').write_text(json.dumps(payload,indent=2))
    print('correlation',k,'partitions',len(partitions),'ordered',payload['ordered_count'],flush=True)
