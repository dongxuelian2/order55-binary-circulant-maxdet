"""Generate exact necessary folded norm thresholds before binary lifting."""
from fractions import Fraction
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
records = [json.loads(x) for x in (ROOT/'results/order55_incumbents.jsonl').read_text().splitlines()]
M = max(max(int(x['determinant']), int(x['complement_determinant'])) for x in records)
rows = []
for m in (5, 11):
    pairs = (m-1)//2
    rest = 27-pairs
    for k in range(24, 28):
        for squares in range(k*k//m, (55//m)*k+1):
            sector_sum = Fraction(m*squares-k*k, 2)
            other_sum = Fraction(k*(55-k), 2)-sector_sum
            if sector_sum <= 0 or other_sum <= 0:
                continue
            sector_max = (sector_sum/pairs)**pairs
            other_max = (other_sum/rest)**rest
            if (55-k)*sector_max*other_max < M:
                continue
            need = Fraction(M, 55-k)/other_max
            lower = -(-need.numerator//need.denominator)
            rows.append((m, k, squares, lower))
(ROOT/'results/order55_profile_thresholds.txt').write_text('\n'.join(' '.join(map(str,row)) for row in rows)+'\n')
(ROOT/'results/order55_profile_thresholds_metadata.json').write_text(json.dumps(dict(incumbent=str(M), rule='(55-k) N_m (S_remaining / remaining_pairs)^remaining_pairs >= M', rows=len(rows)), indent=2))
print('incumbent', M, 'threshold rows', len(rows))
