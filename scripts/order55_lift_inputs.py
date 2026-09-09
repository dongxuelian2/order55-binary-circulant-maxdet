import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
d=json.loads((ROOT/'results/order55_lift_config_k24.json').read_text())
for i,r in enumerate(d['mod5']):
    for j,s in enumerate(d['mod11']):
        values=list(r)+list(s)+d['correlations']
        (ROOT/f'results/order55_lift_k24_{i}_{j}.txt').write_text(' '.join(map(str,values))+'\n')
