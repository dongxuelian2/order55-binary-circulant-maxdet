from pathlib import Path
import json,subprocess,hashlib
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[2];OLD=ROOT/"certificates/order55_global";OUT=ROOT/"certificates/order55_sign/binary_companion"
OUT.mkdir(exist_ok=True);exe=ROOT/"native/order55_sign/independent_lift.exe"
def run(i):
 inp=OLD/f"lift_part{i}.txt";audit=OUT/f"part{i}_audit.jsonl";words=OUT/f"part{i}_words.txt"
 subprocess.run([str(exe),str(inp),str(audit),str(words),"6"],check=True)
 expected={l.split()[0] for l in inp.read_text().splitlines()}
 a=[json.loads(l) for l in audit.read_text().splitlines()];b=[json.loads(l) for l in (OLD/f"lift_part{i}_audit.jsonl").read_text().splitlines()]
 assert len(a)==len(expected) and {r["task"] for r in a}==expected
 assert {r["task"]:(r["joined_words"],r["solutions"]) for r in a}=={r["task"]:(r["joined_words"],r["solutions"]) for r in b}
 assert set(words.read_text().splitlines())==set((OLD/f"lift_part{i}_words.txt").read_text().splitlines())
 return dict(part=i,tasks=len(a),joined_words=sum(r["joined_words"] for r in a),solutions=sum(r["solutions"] for r in a),input_sha256=hashlib.sha256(inp.read_bytes()).hexdigest())
with ThreadPoolExecutor(max_workers=4) as pool:rows=list(pool.map(run,range(8)))
(OUT/"full_fiber_report.json").write_text(json.dumps(dict(status="PASS",scope="Completes the previously missing independent full binary-55 fibers, after sign theorem closure",implementation="audit/order55/independent_torus_lift.cpp",source_sha256=hashlib.sha256((ROOT/"audit/order55/independent_torus_lift.cpp").read_bytes()).hexdigest(),tasks=sum(r["tasks"] for r in rows),joined_words=sum(r["joined_words"] for r in rows),solutions=sum(r["solutions"] for r in rows),parts=rows),indent=2)+"\n")
print("BINARY COMPANION INDEPENDENT FULL FIBERS PASS",sum(r["tasks"] for r in rows),flush=True)

