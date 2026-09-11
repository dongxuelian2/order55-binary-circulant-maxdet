from pathlib import Path
import sys,json
ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/"src"))
from maxdet.sign55 import certify,weight_screen
rows=[json.loads(l) for l in (ROOT/"research/order55_sign/search.jsonl").read_text().splitlines() if l]
best={}
for r in rows:
    if r["k"] not in best or r["score"]>best[r["k"]]["score"]:best[r["k"]]=r
out=ROOT/"certificates/order55_sign";out.mkdir(exist_ok=True,parents=True)
certificates=[]
for k,r in sorted(best.items()):
    c=certify(r["word"]);certificates.append(c)
    print(k,c["normalized"],c["word"],flush=True)
(out/"incumbents.json").write_text(json.dumps(certificates,indent=2)+"\n")
winner=max(certificates,key=lambda c:int(c["normalized"]))
(out/"incumbent.json").write_text(json.dumps(winner,indent=2)+"\n")
screen=weight_screen(int(winner["normalized"]))
(out/"weight_screen.json").write_text(json.dumps(screen,indent=2)+"\n")
print("BEST",winner["normalized"],"weight",winner["weight"])
print("ACTIVE",[(r["k"],r.get("maximum_half_squares"),r["upper_floor"]) for r in screen if not r["excluded"]])

