"""Prepare exact margin orientations for a particular correlation target."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
k=24
target=json.loads((ROOT/f'results/order55_correlation_targets_k{k}.jsonl').read_text(encoding='utf-8-sig').splitlines()[0])
c=[k]+target['correlations']+list(reversed(target['correlations']))
out=dict(k=k,correlations=c)
for m in (5,11):
    want=tuple(sum(c[i::m]) for i in range(m))
    good=set()
    for line in (ROOT/f'results/order55_profiles_m{m}_k{k}.jsonl').read_text(encoding='utf-8-sig').splitlines():
        a=json.loads(line)['profile']
        for u in range(1,m):
            for t in range(m):
                b=tuple(a[(u*i+t)%m] for i in range(m))
                h=tuple(sum(b[i]*b[(i+s)%m] for i in range(m)) for s in range(m))
                if h==want:good.add(b)
    # Translation preserves c: one lexicographically minimal rotation suffices
    # independently in both CRT coordinates, but unit orientations stay fixed.
    good={min(b[t:]+b[:t] for t in range(m)) for b in good}
    out[f'mod{m}']=sorted(good)
    print(m,len(good),'oriented margin necklaces',sorted(good),flush=True)
(ROOT/'results/order55_lift_config_k24.json').write_text(json.dumps(out,indent=2))
# Torus coordinate (i,j) corresponds to t=11*i+45*j mod 55.
out['torus_target']=[[c[(11*i+45*j)%55] for j in range(11)] for i in range(5)]
(ROOT/'results/order55_lift_config_k24.json').write_text(json.dumps(out,indent=2))
