"""Prepare NEW sign threshold profile universe. Does not assert completeness of lifts."""
from pathlib import Path
import sys,json,hashlib,subprocess
from fractions import Fraction as F
from collections import defaultdict
from functools import lru_cache
from concurrent.futures import ThreadPoolExecutor
ROOT=Path(__file__).resolve().parents[2]
sys.path[:0]=[str(ROOT/"src"),str(ROOT/"scripts")]
from maxdet.sign55 import weight_screen,fourth_upper
from maxdet.order55 import balanced_squares,PRIMES,moment_product_upper
from maxdet.boxed_moment import boxed
from build_order55_global import partitions,cosine_table
OUT=ROOT/"certificates/order55_sign"
M=int(json.loads((OUT/"incumbent.json").read_text())["normalized"])
def write(name,data):(OUT/name).write_text(json.dumps(data,indent=2)+"\n")
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
weights=weight_screen(M);write("weight_screen.json",weights)
active=[r["k"] for r in weights if not r["excluded"]]
write("global.json",dict(status="INCOMPLETE",incumbent=str(M),active_weights=active,global_completeness="NOT PROVED",cleanroom_profile_audit="PENDING",cleanroom_lift_audit="PENDING"))
cos=cosine_table()
(OUT/"cosine_intervals.txt").write_text("\n".join(f"{a} {b}" for a,b in cos)+"\n")
thresholds=[]
for m in (5,11):
    h=(m-1)//2;rest=27-h;cap=55//m
    assert F(m*m*cap*cap,m-1)**h<PRIMES[0]
    for k in active:
        for square in range(k*k//m,cap*k+1):
            T=F(m*square-k*k,2);R=F(k*(55-k),2)-T
            if T<=0 or R<=0 or (55-2*k)*(T/h)**h*(R/rest)**rest<M:continue
            need=F(M,55-2*k)/(R/rest)**rest
            thresholds.append((m,k,square,-(-need.numerator//need.denominator)))
(OUT/"profile_thresholds.txt").write_text("\n".join(" ".join(map(str,r)) for r in thresholds)+"\n")
# Fold and lift programs do not reference the objective.
for name in ("order55_profiles","order55_torus_lift"):
    dest=ROOT/f"native/order55_sign/{name}.exe"
    subprocess.run(["clang++","-O3","-std=c++20","-march=native",str(ROOT/f"native/{name}.cpp"),"-o",str(dest)],check=True)
def folds(task):
    k,m=task;path=OUT/f"profiles_m{m}_k{k}.jsonl"
    with path.open("w") as out:
        run=subprocess.run([str(ROOT/"native/order55_sign/order55_profiles.exe"),str(OUT/"profile_thresholds.txt"),str(m),str(k)],stdout=out,stderr=subprocess.PIPE,text=True,check=True)
    rows=[json.loads(l) for l in path.read_text().splitlines()]
    lookup=defaultdict(set)
    for row in rows:
        a=row["profile"]
        for u in range(1,m):
            b=tuple(a[u*i%m] for i in range(m));b=min(b[t:]+b[:t] for t in range(m))
            h=tuple(sum(b[t]*b[(t+s)%m] for t in range(m)) for s in range(m))
            lookup[h].add(b)
    sigs=sorted({h[:(m+1)//2] for h in lookup})
    (OUT/f"foldcorr_m{m}_k{k}.txt").write_text("\n".join(" ".join(map(str,s)) for s in sigs)+"\n")
    print("FOLDS",k,m,len(rows),"signatures",len(sigs),flush=True)
    return dict(k=k,m=m,survivors=len(rows),signatures=len(sigs),sha256=digest(path),stderr=run.stderr)
with ThreadPoolExecutor(max_workers=6) as pool:audits=list(pool.map(folds,[(k,m) for k in active for m in (5,11)]))
write("fold_manifest.json",audits)
moment=lru_cache(None)(moment_product_upper)
fourth_upper=lru_cache(None)(fourth_upper)
# Rigorous joint-fold bounds. These use every NEW fold survivor.
joint=[]
for k in active:
    f={m:[json.loads(l) for l in (OUT/f"profiles_m{m}_k{k}.jsonl").read_text().splitlines()] for m in (5,11)}
    survivors=[];best=0
    for r in f[5]:
        for s in f[11]:
            A=(r["squares"]-k)//2;B=(s["squares"]-k)//2;C=k*(k-1)//2-A-B
            if min(A,B,C)<0:continue
            lower=balanced_squares(A,5)+balanced_squares(B,2)+balanced_squares(C,20)
            T=F(55*k+k*k-5*r["squares"]-11*s["squares"],2)
            if T<=0:continue
            # First sector bound plus global fourth moment.
            u=min((55-2*k)*r["norm"]*s["norm"]*(T/20)**20,fourth_upper(k,lower))
            best=max(best,int(u))
            if u>=M:survivors.append(dict(mod5=r["profile"],mod11=s["profile"],N5=r["norm"],N11=s["norm"],upper_floor=int(u),minimum_half_squares=lower))
    write(f"joint_k{k}.json",survivors)
    joint.append(dict(k=k,pairs=len(f[5])*len(f[11]),survivors=len(survivors),upper_floor=best,excluded=not survivors))
    print("JOINT",joint[-1],flush=True)
write("joint_manifest.json",joint)
manifest=[]
for k in active:
    if next(r for r in joint if r["k"]==k)["excluded"]:continue
    limit=weights[k-1]["maximum_half_squares"]
    for idx,p in enumerate(partitions(k,limit)):
        base=k*(k-1)//54;ds=sorted(v-base for v in p["values"]);caps=[]
        for j in (1,5,11):
            cs=sorted(cos[min(j*s%55,55-j*s%55)] for s in range(1,28))
            caps.append(F(k-base)+sum(F(v*(pair[1] if v>=0 else pair[0]),1<<24) for v,pair in zip(ds,cs)))
        cap=max(caps);S1=F(k*(55-k),2);S2=F(55*(k*k+2*p["squares"])-k**4,2)
        upper=int((55-2*k)*boxed(S1,S2,27,cap))
        assert (2*(55-2*k))**2*S2**27 < (PRIMES[0]*PRIMES[1])**2*27**27
        p.update(k=k,partition=idx,cap=str(cap),upper_floor=upper,excluded=upper<M)
        manifest.append(p)
        if upper>=M:
            (OUT/f"corr_k{k}_part{idx}.txt").write_text(" ".join(map(str,p["values"]))+"\n")
    relevant=[p for p in manifest if p["k"]==k and not p["excluded"]]
    print("PARTITIONS",k,len(relevant),"formal profiles",sum(p["ordered_count"] for p in relevant),flush=True)
write("profile_manifest.json",dict(status="PREPARED; ordered enumeration pending",incumbent=str(M),partitions=manifest))


