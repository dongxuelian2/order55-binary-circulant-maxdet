"""Run every sign-specific retained correlation multiset, then materialize all lifts."""
from pathlib import Path
import sys,json,hashlib,subprocess,time
from fractions import Fraction as F
from math import isqrt
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
OUT=ROOT/"certificates/order55_sign"
def write(name,data):(OUT/name).write_text(json.dumps(data,indent=2)+"\n")
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
manifest=json.loads((OUT/"profile_manifest.json").read_text())
M=int(manifest["incumbent"])
parts=[p for p in manifest["partitions"] if not p["excluded"]]
cos=[tuple(map(int,l.split())) for l in (OUT/"cosine_intervals.txt").read_text().splitlines()]
for p in parts:
 k=p["k"];base=k*(k-1)//54;ds=[v-base for v in p["values"]]
 S2=F(55*(k*k+2*p["squares"])-k**4,2)
 err=F(27*sum(abs(v) for v in ds)*max(b-a for a,b in cos),1<<24)
 A=F(isqrt(int(27*S2))+1)+err
 prefix=max([F(1)]+[(A/(16*h))**h for h in range(1,28)])
 assert ((1<<64)+27)*prefix*A*(1<<24)+(1<<28)<1<<128
source=(ROOT/"native/order55_correlation_screen.cpp").read_text()
assert source.count("55-k")==2
source=source.replace("55-k","55-2*k")
native=ROOT/"native/order55_sign/correlation_screen.cpp"
native.write_text("// Sign objective adaptation: normalized determinant coefficient n-2k.\n"+source)
exe=native.with_suffix(".exe")
subprocess.run(["clang++","-O3","-std=c++20","-march=native","-Wall","-Wextra",str(native),"-o",str(exe)],check=True)
def run(p):
 k,idx=p["k"],p["partition"];stem=f"corr_k{k}_part{idx}"
 path=OUT/f"{stem}_targets.jsonl";auditpath=OUT/f"{stem}_audit.json"
 cmd=[str(exe),str(k),str(M),str(OUT/f"{stem}.txt"),str(OUT/"foldcorr_m"),str(auditpath),str(OUT/"cosine_intervals.txt")]
 with path.open("w") as f:subprocess.run(cmd,stdout=f,check=True)
 audit=json.loads(auditpath.read_text())
 assert audit["ordered_profiles"]==p["ordered_count"]
 print("PROFILE COMPLETE",k,idx,audit,flush=True)
 return dict(k=k,partition=idx,audit=audit,target_file=path.name,target_sha256=sha(path),source_sha256=sha(native))
audits=[]
with ThreadPoolExecutor(max_workers=6) as pool:
 for future in as_completed([pool.submit(run,p) for p in parts]):
  audits.append(future.result());write("profile_progress.json",audits)
audits.sort(key=lambda r:(r["k"],r["partition"]))
targets=[]
for r in audits:
 for line in (OUT/r["target_file"]).read_text().splitlines():
  t=json.loads(line);t["source_partition"]=r["partition"];targets.append(t)
targets.sort(key=lambda t:(-int(t["absolute_profile_product"]),t["k"],t["correlations"]))
def target_hash(rows):
 return hashlib.sha256("".join(f'{r["k"]} '+",".join(map(str,r["correlations"]))+"\n" for r in sorted(rows,key=lambda r:(r["k"],r["correlations"]))).encode()).hexdigest()
write("profile_coverage.json",dict(status="PRODUCTION PROFILE GENERATION COMPLETE; independent audit pending",incumbent=str(M),total_formal_profiles=sum(p["ordered_count"] for p in manifest["partitions"]),explicitly_enumerated=sum(p["ordered_count"] for p in parts),boxed_pruned=sum(p["ordered_count"] for p in manifest["partitions"] if p["excluded"]),surviving_profiles=len(targets),target_set_sha256=target_hash(targets),audits=audits))
write("correlation_targets.json",targets)
lookups={}
for k in sorted({t["k"] for t in targets}):
 for m in (5,11):
  lookup=defaultdict(set)
  for line in (OUT/f"profiles_m{m}_k{k}.jsonl").read_text().splitlines():
   a=json.loads(line)["profile"]
   for u in range(1,m):
    b=tuple(a[u*i%m] for i in range(m));b=min(b[t:]+b[:t] for t in range(m))
    h=tuple(sum(b[t]*b[(t+s)%m] for t in range(m)) for s in range(m))
    lookup[h].add(b)
  lookups[k,m]=lookup
workers=8;handles=[(OUT/f"lift_part{i}.txt").open("w") for i in range(workers)]
coverage=[];total=0
for idx,t in enumerate(targets):
 k=t["k"];c=[k]+t["correlations"]+list(reversed(t["correlations"]))
 rs=sorted(lookups[k,5][tuple(sum(c[i::5]) for i in range(5))]);ss=sorted(lookups[k,11][tuple(sum(c[i::11]) for i in range(11))])
 for i,r in enumerate(rs):
  for j,s in enumerate(ss):
   task=f"{idx}_{i}_{j}";handles[total%workers].write(task+" "+" ".join(map(str,list(r)+list(s)+c))+"\n");total+=1
 coverage.append(dict(id=idx,target=t,row_profiles=rs,column_profiles=ss,tasks=len(rs)*len(ss)))
for h in handles:h.close()
write("lift_targets.json",coverage)
write("lift_manifest.json",dict(status="INPUTS COMPLETE; binary enumeration pending",targets=len(targets),tasks=total,partitions={f"lift_part{i}.txt":sha(OUT/f"lift_part{i}.txt") for i in range(workers)}))
print("ALL SIGN PROFILES",len(targets),"LIFT TASKS",total,flush=True)

