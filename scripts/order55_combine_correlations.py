from pathlib import Path
import json
root=Path('results')
boxed={(r['k'],r['partition']):r for r in json.loads((root/'order55_boxed_moment_screen.json').read_text())}
for k in (26,27):
    targets=[];coverage=[];complete=True
    parts=json.loads((root/f'order55_correlation_partitions_k{k}.json').read_text())['partitions']
    for idx,p in enumerate(parts):
        b=boxed[k,idx];audit=root/f'order55_correlation_work/k{k}_part{idx:02}_audit.json'
        if b['excluded']:
            coverage.append(dict(partition=idx,route='boxed moment',upper=b['upper'],ordered_count=p['ordered_count']))
        elif audit.exists():
            a=json.loads(audit.read_text());assert a['ordered_profiles']==p['ordered_count']
            coverage.append(dict(partition=idx,route='exact enumeration',audit=str(audit),ordered_count=p['ordered_count']))
            target=root/f'order55_correlation_work/k{k}_part{idx:02}_targets.jsonl'
            targets.extend(json.loads(x) for x in target.read_text().splitlines())
        else:
            complete=False;coverage.append(dict(partition=idx,route='PENDING'));print('PENDING',k,idx)
    (root/f'order55_correlation_combined_k{k}.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in targets))
    (root/f'order55_correlation_coverage_k{k}.json').write_text(json.dumps(dict(complete=complete,targets=len(targets),coverage=coverage),indent=2))
    print(k,'complete',complete,'targets',len(targets))
