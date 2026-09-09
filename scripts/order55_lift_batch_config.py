"""Complete margin tasks for every requested correlation target, hash bound."""
from collections import defaultdict
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]
k=int(sys.argv[1]);source=Path(sys.argv[2]);prefix=ROOT/f'results/order55_lift_batch_k{k}'
lookups={}
for m in (5,11):
    lookup=defaultdict(set)
    for line in (ROOT/f'results/order55_profiles_m{m}_k{k}.jsonl').read_text(encoding='utf-8-sig').splitlines():
        a=json.loads(line)['profile']
        for u in range(1,m):
            b=tuple(a[u*i%m] for i in range(m));b=min(b[t:]+b[:t] for t in range(m))
            h=tuple(sum(b[t]*b[(t+s)%m] for t in range(m)) for s in range(m))
            lookup[h].add(b)
    lookups[m]=lookup
rows=[];total=0;files=[Path(str(prefix)+f'_part{i}.txt').open('w') for i in range(8)]
targets=[json.loads(x) for x in source.read_text(encoding='utf-8-sig').splitlines()]
targets.sort(key=lambda x:int(x['absolute_profile_product']),reverse=True)
for idx,t in enumerate(targets):
    c=[k]+t['correlations']+list(reversed(t['correlations']))
    rs=sorted(lookups[5][tuple(sum(c[i::5]) for i in range(5))])
    ss=sorted(lookups[11][tuple(sum(c[i::11]) for i in range(11))])
    for i,r in enumerate(rs):
        for j,s in enumerate(ss):
            task=f'{idx}_{i}_{j}'
            files[total%8].write(task+' '+' '.join(map(str,list(r)+list(s)+c))+'\n');total+=1
    rows.append(dict(id=idx,target=t,row_profiles=len(rs),column_profiles=len(ss),tasks=len(rs)*len(ss)))
for f in files:f.close()
meta=dict(k=k,source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),targets=len(targets),tasks=total,coverage=rows)
Path(str(prefix)+'_manifest.json').write_text(json.dumps(meta))
print('targets',len(targets),'complete margin tasks',total,flush=True)
