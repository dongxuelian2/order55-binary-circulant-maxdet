from pathlib import Path
import json
root=Path('results')
for k in (25,26,27):
    manifest=root/f'order55_lift_batch_k{k}_manifest.json'
    if not manifest.exists():continue
    d=json.loads(manifest.read_text());seen=set();joined=0;solutions=0
    for p in root.glob(f'order55_lift_batch_k{k}_part*_audit.jsonl'):
        for line in p.read_text().splitlines():
            try:a=json.loads(line)
            except json.JSONDecodeError:continue
            assert a['task'] not in seen;seen.add(a['task']);joined+=a['joined_words'];solutions+=a['solutions']
    print(k,'completed',len(seen),'of',d['tasks'],'joined',joined,'solutions',solutions)
