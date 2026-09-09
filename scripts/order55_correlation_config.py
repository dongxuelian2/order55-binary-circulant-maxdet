import json
from pathlib import Path
from itertools import product
from maxdet.order55 import PRIMES
ROOT=Path(__file__).resolve().parents[1]
for k in range(24,28):
    parts=json.loads((ROOT/f'results/order55_correlation_partitions_k{k}.json').read_text())
    (ROOT/f'results/order55_correlation_parts_k{k}.txt').write_text('\n'.join(' '.join(map(str,p['values'])) for p in parts['partitions'])+'\n')
    # All real spectral values (even for non-PSD profiles) have this fixed
    # sum of squares. AM-GM applied to q_j^2 certifies absolute CRT capacity.
    max_s2=(55*(k*k+2*parts['maximum_half_squares'])-k**4)//2
    assert (2*(55-k))**2 * max_s2**27 < (PRIMES[0]*PRIMES[1])**2 * 27**27
    for m in (5,11):
        profiles=[json.loads(x)['profile'] for x in (ROOT/f'results/order55_profiles_m{m}_k{k}.jsonl').read_text(encoding='utf-8-sig').splitlines()]
        signatures=set()
        for a in profiles:
            h=[sum(a[t]*a[(t+s)%m] for t in range(m)) for s in range(m)]
            for u in range(1,m):signatures.add(tuple(h[u*s%m] for s in range((m+1)//2)))
        (ROOT/f'results/order55_foldcorr_m{m}_k{k}.txt').write_text('\n'.join(' '.join(map(str,s)) for s in sorted(signatures))+'\n')
        print(k,m,'folded correlation signatures',len(signatures))
