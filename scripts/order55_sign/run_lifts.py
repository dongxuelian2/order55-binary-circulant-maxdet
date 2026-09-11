from pathlib import Path
import sys,json,subprocess,time,hashlib
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from maxdet.sign55 import certify
OUT=ROOT/"certificates/order55_sign"
def write(n,d):(OUT/n).write_text(json.dumps(d,indent=2)+"\n")
mode=sys.argv[1] if len(sys.argv)>1 else "production"
independent=mode=="independent";prefix="independent" if independent else "production"
source=ROOT/("audit/order55/independent_torus_lift.cpp" if independent else "native/order55_torus_lift.cpp")
exe=ROOT/f"native/order55_sign/{prefix}_lift.exe"
subprocess.run(["clang++","-O3","-std=c++20","-march=native","-Wall","-Wextra",str(source),"-o",str(exe)],check=True)
def run(i):
 inp=OUT/f"lift_part{i}.txt";audit=OUT/f"{prefix}_lift_part{i}_audit.jsonl";words=OUT/f"{prefix}_lift_part{i}_words.txt"
 start=time.time()
 if independent:
  cmd=[str(exe),str(inp),str(audit),str(words),"6"];subprocess.run(cmd,check=True)
 else:
  cmd=[str(exe),str(inp),"5",str(audit),"batch"]
  with words.open("w") as out:subprocess.run(cmd,stdout=out,check=True)
 expected={l.split()[0] for l in inp.read_text().splitlines()}
 rows=[json.loads(l) for l in audit.read_text().splitlines()]
 assert len(rows)==len(expected) and {r["task"] for r in rows}==expected
 report=dict(part=i,tasks=len(rows),joined_words=sum(r["joined_words"] for r in rows),solutions=sum(r["solutions"] for r in rows),seconds=time.time()-start,source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),audit_sha256=hashlib.sha256(audit.read_bytes()).hexdigest(),words_sha256=hashlib.sha256(words.read_bytes()).hexdigest())
 print(prefix,report,flush=True);return report
reports=[]
with ThreadPoolExecutor(max_workers=8 if not independent else 4) as pool:
 for f in as_completed([pool.submit(run,i) for i in range(8)]):
  reports.append(f.result());write(f"{prefix}_lift_progress.json",reports)
write(f"{prefix}_lift_report.json",dict(status="FULL LISTED FIBERS ENUMERATED",parts=sorted(reports,key=lambda r:r["part"]),total_tasks=sum(r["tasks"] for r in reports),joined_words=sum(r["joined_words"] for r in reports),solutions=sum(r["solutions"] for r in reports),split=6 if independent else 5))
if independent:
 for i in range(8):
  if not (OUT/f"production_lift_part{i}_words.txt").exists():continue
  a=set((OUT/f"production_lift_part{i}_words.txt").read_text().splitlines());b=set((OUT/f"independent_lift_part{i}_words.txt").read_text().splitlines())
  assert a==b
else:
 words=[]
 for i in range(8):
  for l in (OUT/f"production_lift_part{i}_words.txt").read_text().splitlines():
   task,w=l.split();words.append(dict(task=task,word=w))
 write("lifted_words.json",words)
 certs={}
 for w in words:
  c=certify(w["word"]);certs[c["word"]]=c
 assert certs,"incumbent must lift"
 best=max(int(c["normalized"]) for c in certs.values())
 write("lifted_exact_candidates.json",list(certs.values()))
 write("candidate_winners.json",[c for c in certs.values() if int(c["normalized"])==best])
 print("EXACT LIFT BEST",best,"CLASSES",sum(int(c["normalized"])==best for c in certs.values()),flush=True)

