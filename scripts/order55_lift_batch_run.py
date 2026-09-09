from pathlib import Path
from concurrent.futures import ThreadPoolExecutor,as_completed
import subprocess,sys
ROOT=Path(__file__).resolve().parents[1]
k=int(sys.argv[1]);prefix=ROOT/f'results/order55_lift_batch_k{k}'
def run(i):
    cmd=[str(ROOT/'native/order55_torus_lift.exe'),str(prefix)+f'_part{i}.txt','5',str(prefix)+f'_part{i}_audit.jsonl','batch']
    with open(str(prefix)+f'_part{i}_words.txt','w') as f:subprocess.run(cmd,stdout=f,check=True)
    print('COMPLETE',k,i,flush=True)
with ThreadPoolExecutor(max_workers=4) as pool:
    for f in as_completed([pool.submit(run,i) for i in range(8)]):f.result()
print('ALL LIFT TASKS COMPLETE',k,flush=True)
