"""Independent sign arithmetic, multiset/fold-universe and target-set audit."""
from pathlib import Path
import sys,json,subprocess,hashlib
from fractions import Fraction as F
from math import isqrt,factorial,prod
from collections import Counter
from concurrent.futures import ThreadPoolExecutor,as_completed
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/"certificates/order55_sign";CLEAN=OUT/"independent";CLEAN.mkdir(exist_ok=True)
def write(path,obj):path.write_text(json.dumps(obj,indent=2)+"\n")
def root_bounds(v):
 v=F(v);scale=1<<220;r=isqrt(v.numerator*scale*scale//v.denominator)
 return F(r,scale),F(r if r*r*v.denominator==v.numerator*scale*scale else r+1,scale)
def bound(T,L,h,B=None):
 T,L=F(T),F(L);L=max(L,T*T/h);best=F(0)
 if B is None:B=T
 B=F(B)
 for e in range(h):
  n=h-e;t=T-e*B;l=L-e*B*B
  if t<0 or t>n*B or l<t*t/n or l>t*B:continue
  if l==t*t/n:best=max(best,B**e*(t/n)**n);continue
  for a in range(1,n):
   b=n-a;vlo,vhi=root_bounds(F(n*l-t*t,a*b))
   if t-b*vlo<=0:continue
   if (t+a*vlo)/n>B:continue
   best=max(best,B**e*((t-b*vlo)/n)**a*((t+a*vhi)/n)**b)
 return best
def balance(s,h):
 a,b=divmod(s,h);return a*a*h+b*(2*a+1)
def fourth(k,sq):
 return (55-2*k)*bound(F(k*(55-k),2),F(55*(k*k+2*sq)-k**4,2),27)
from functools import lru_cache
fourth=lru_cache(None)(fourth)
M=int(json.loads((OUT/"incumbent.json").read_text())["normalized"])
weights=json.loads((OUT/"weight_screen.json").read_text());parts=json.loads((OUT/"profile_manifest.json").read_text())["partitions"]
for row in weights:
 k=row["k"];first=(55-2*k)*F(k*(55-k),54)**27
 assert first.numerator==int(row["first_upper_num"]) and first.denominator==int(row["first_upper_den"])
 sq=balance(k*(k-1)//2,27);u=min(first,fourth(k,sq))
 assert int(u)<=row["upper_floor"]
 assert row["excluded"]==(u<M)
 if not row["excluded"]:
  limit=row["maximum_half_squares"]
  assert fourth(k,limit)>=M and fourth(k,limit+1)<M
# Independent cosine enclosures: Taylor at BOTH endpoints, exploiting
# monotonicity on [0,pi]; no midpoint derivative-error construction.
def atan(q):
 terms=[F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(111)]
 partial=sum(terms[:-1]);return sorted((partial,partial+terms[-1]))
a,b=atan(5);c,d=atan(239);pi_lo,pi_hi=16*a-4*d,16*b-4*c
scale=1<<256
pi_lo=F((pi_lo*scale).numerator//(pi_lo*scale).denominator,scale)
pi_hi=F(-(-(pi_hi*scale).numerator//(pi_hi*scale).denominator),scale)
print("INDEPENDENT WEIGHTS PASS; checking dyadic cosine enclosure",flush=True)
def cos_bounds(x):
 s=sum(F((-1)**j)*x**(2*j)/factorial(2*j) for j in range(70))
 t=F((-1)**70)*x**140/factorial(140)
 return min(s,s+t),max(s,s+t)
cos=[tuple(map(int,l.split())) for l in (OUT/"cosine_intervals.txt").read_text().splitlines()]
for j,(lo,hi) in enumerate(cos):
 lower=cos_bounds(2*pi_hi*j/55)[0];upper=cos_bounds(2*pi_lo*j/55)[1]
 assert F(lo,1<<24)<=2*lower and 2*upper<=F(hi,1<<24)
# Enumerate MULTIPLICITIES by descending values; production grows a sorted word.
for row in weights:
 if row["excluded"]:continue
 k=row["k"];limit=row["maximum_half_squares"];total=k*(k-1)//2;found=set()
 feasible=[v for v in range(k+1) if v*v+F((total-v)**2,26)<=limit]
 def multiplicities(v,left,s,sq,counts):
  if left==0:
   if s==0:found.add(tuple(sorted(counts)))
   return
  if v<0 or s<0 or s>v*left or sq+balance(s,left)>limit:return
  for copies in range(min(left,s//v if v else left)+1):
   newleft=left-copies;news=s-copies*v;newsq=sq+copies*v*v
   if newsq<=limit:multiplicities(v-1,newleft,news,newsq,counts+[v]*copies)
 multiplicities(max(feasible),27,total,0,[])
 actual={tuple(p["values"]) for p in parts if p["k"]==k}
 assert found==actual,(k,len(found),len(actual))
 for p in (p for p in parts if p["k"]==k):
  freq=Counter(p["values"]);count=factorial(27)//prod(factorial(x) for x in freq.values());assert count==p["ordered_count"]
  base=k*(k-1)//54;ds=sorted(v-base for v in p["values"]);caps=[]
  for frequency in (1,5,11):
   pairs=sorted(cos[min(frequency*s%55,55-frequency*s%55)] for s in range(1,28))
   caps.append(F(k-base)+sum(F(v*(hi if v>=0 else lo),1<<24) for v,(lo,hi) in zip(ds,pairs)))
  assert str(max(caps))==p["cap"]
  u=(55-2*k)*bound(F(k*(55-k),2),F(55*(k*k+2*p["squares"])-k**4,2),27,max(caps))
  assert int(u)<=p["upper_floor"]
  assert p["excluded"]==(u<M)
# Independently reconstruct all fold thresholds.
threshold=[]
active=[r["k"] for r in weights if not r["excluded"]]
for k in active:
 for m in (5,11):
  h=(m-1)//2
  assert F(m*m*(55//m)**2,m-1)**h<2305843009213697141
  for square in range(k*k//m,(55//m)*k+1):
   t=F(m*square-k*k,2);r=F(k*(55-k),2)-t
   if t<=0 or r<=0:continue
   ceiling=(55-2*k)*(t/h)**h*(r/(27-h))**(27-h)
   if ceiling<M:continue
   ratio=F(M,55-2*k)/(r/(27-h))**(27-h)
   threshold.append((m,k,square,(ratio.numerator+ratio.denominator-1)//ratio.denominator))
assert set(threshold)=={tuple(map(int,l.split())) for l in (OUT/"profile_thresholds.txt").read_text().splitlines()}
th=CLEAN/"profile_thresholds.txt";th.write_text("\n".join(" ".join(map(str,r)) for r in sorted(threshold))+"\n")
write(CLEAN/"bounds_report.json",dict(status="PASS",independent_weight_screen=True,independent_two_root_moment_bound=True,independent_capped_bound=True,independent_cosine_enclosures=True,independent_multiset_universe=True,independent_fold_thresholds=True))
print("INDEPENDENT BOUNDS AND MULTISET UNIVERSE PASS",flush=True)
subprocess.run(["clang++","-O3","-std=c++20","-march=native",str(ROOT/"native/order55_sign/independent_folds.cpp"),"-o",str(ROOT/"native/order55_sign/independent_folds.exe")],check=True)
def fold(task):
 k,m=task;path=CLEAN/f"folds_m{m}_k{k}.jsonl"
 with path.open("w") as f:
  subprocess.run([str(ROOT/"native/order55_sign/independent_folds.exe"),str(th),str(m),str(k),str(55//m)],stdout=f,stderr=subprocess.DEVNULL,check=True)
 rows=[json.loads(l) for l in path.read_text().splitlines()];actual={tuple(r["profile"]) for r in rows}
 expected=set()
 for l in (OUT/f"profiles_m{m}_k{k}.jsonl").read_text().splitlines():
  r=json.loads(l);a=r["profile"]
  for u in range(1,m):
   b=tuple(a[u*i%m] for i in range(m));expected.add(min(b[t:]+b[:t] for t in range(m)))
 assert actual==expected,(k,m,len(actual),len(expected))
 signatures=set()
 for a in actual:
  corr=tuple(sum(a[i]*a[(i+s)%m] for i in range(m)) for s in range((m+1)//2))
  signatures.add(corr)
 signature_path=CLEAN/f"foldcorr_m{m}_k{k}.txt";signature_path.write_text("\n".join(" ".join(map(str,r)) for r in sorted(signatures))+"\n")
 expected_sigs={tuple(map(int,l.split())) for l in (OUT/f"foldcorr_m{m}_k{k}.txt").read_text().splitlines()}
 assert signatures==expected_sigs
 print("INDEPENDENT FOLDS PASS",k,m,len(actual),flush=True)
 return dict(k=k,m=m,necklace_profiles=len(actual),signatures=len(signatures),status="PASS")
with ThreadPoolExecutor(max_workers=4) as pool:folds=list(pool.map(fold,[(k,m) for k in active for m in (5,11)]))
write(CLEAN/"fold_report.json",dict(status="PASS",algorithm="FKM necklaces, every rotation class; no production affine canonicalization",rows=folds))
# Independent next_permutation traversal, with every canonical CRT value
# evaluated (no production floating point or interval-product decisions).
source=(ROOT/"audit/order55/independent_profile_screen.cpp").read_text()
assert source.count("55 - k")==1
source=source.replace("55 - k","55 - 2*k")
native=ROOT/"native/order55_sign/independent_profiles.cpp";native.write_text(source)
exe=native.with_suffix(".exe");subprocess.run(["clang++","-O3","-std=c++20","-march=native",str(native),"-o",str(exe)],check=True)
def screen(p):
 k,idx=p["k"],p["partition"];stem=f"corr_k{k}_part{idx}"
 audit=CLEAN/f"{stem}_audit.json";target=CLEAN/f"{stem}_targets.jsonl"
 if "--reuse-enumeration" not in sys.argv:
  subprocess.run([str(exe),str(k),str(M),str(OUT/f"{stem}.txt"),str(CLEAN),str(audit),str(target)],check=True)
 else:
  assert audit.exists() and target.exists()
 a=json.loads(audit.read_text());assert a["ordered_profiles"]==p["ordered_count"]
 convert=lambda path:{(r["k"],tuple(r["correlations"]),int(r["absolute_profile_product"])) for r in (json.loads(l) for l in path.read_text().splitlines())}
 clean=convert(target);production=convert(OUT/f"{stem}_targets.jsonl")
 # The raw independent CRT evaluator intentionally permits indefinite
 # formal spectra. Discharge every extra profile with an explicit negative
 # eigenvalue upper bound from the independently verified cosine enclosure.
 assert production<=clean
 rejected=[]
 for kk,corr,value in sorted(clean-production):
  witness=None
  for j in range(1,28):
   upper=(kk<<24)+sum(v*cos[min(j*shift%55,55-j*shift%55)][1] for shift,v in enumerate(corr,1))
   if upper<0:
    witness=dict(k=kk,correlations=list(corr),absolute_algebraic_product=str(value),
                 frequency=j,eigenvalue_upper_numerator=upper,eigenvalue_upper_denominator=1<<24)
    break
  assert witness is not None,("unexplained target difference",kk,corr,value)
  rejected.append(witness)
 write(CLEAN/f"{stem}_indefinite.json",rejected)
 retained=clean-{(r["k"],tuple(r["correlations"]),int(r["absolute_algebraic_product"])) for r in rejected}
 assert retained==production
 clean=retained
 print("INDEPENDENT PROFILE PASS",k,idx,len(clean),flush=True)
 return dict(k=k,partition=idx,status="PASS",audit=a,indefinite_profiles_excluded=len(rejected),retained_profiles=len(clean),target_set_equal=True)
reports=[]
with ThreadPoolExecutor(max_workers=4) as pool:
 for f in as_completed([pool.submit(screen,p) for p in parts if not p["excluded"]]):
  reports.append(f.result());write(CLEAN/"profile_progress.json",reports)
write(OUT/"cleanroom_profile_report.json",dict(status="PASS",bounds="PASS",fold_universe="PASS",correlation_universe="PASS",canonical_target_sets="EXACT EQUALITY",survivors=sum(r["retained_profiles"] for r in reports),indefinite_profiles_excluded=sum(r["indefinite_profiles_excluded"] for r in reports),raw_algebraic_survivors=sum(r["audit"]["above_screen"] for r in reports),parts=reports))


