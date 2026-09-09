"""Deterministic independent process work units for the correlation screen."""
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib,json,subprocess,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
work=ROOT/'results/order55_correlation_work'
work.mkdir(exist_ok=True)
M='103755181499155107923236157710676'
jobs=[]
for k in (26,27):
    lines=(ROOT/f'results/order55_correlation_parts_k{k}.txt').read_text().splitlines()
    for i,line in enumerate(lines):
        part=work/f'k{k}_part{i:02}.txt'
        part.write_text(line+'\n')
        jobs.append((k,i,part))
def run(job):
    k,i,part=job
    target=work/f'k{k}_part{i:02}_targets.jsonl'
    audit=work/f'k{k}_part{i:02}_audit.json'
    cmd=[str(ROOT/'native/order55_correlation_screen.exe'),str(k),M,str(part),str(ROOT/'results/order55_foldcorr_m'),str(audit),str(ROOT/'results/order55_cosine_intervals.txt')]
    start=time.time()
    with target.open('w') as f: subprocess.run(cmd,stdout=f,check=True,cwd=ROOT)
    row=json.loads(audit.read_text())
    row.update(command=cmd,input_sha256=hashlib.sha256(part.read_bytes()).hexdigest(),output_sha256=hashlib.sha256(target.read_bytes()).hexdigest(),source_sha256=hashlib.sha256((ROOT/'native/order55_correlation_screen.cpp').read_bytes()).hexdigest())
    audit.write_text(json.dumps(row,indent=2))
    print(k,i,row['ordered_profiles'],row['above_screen'],round(time.time()-start,2),flush=True)
    return row
with ThreadPoolExecutor(max_workers=8) as pool:
    futures=[pool.submit(run,job) for job in jobs]
    for f in as_completed(futures):f.result()
print('ALL CORRELATION WORK UNITS COMPLETE',flush=True)
