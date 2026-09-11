from pathlib import Path
import sys,json,subprocess,hashlib
from fractions import Fraction as F
from collections import Counter
from math import factorial,isqrt
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from maxdet.sign55 import first_upper,fourth_upper,bareiss,circulant
from maxdet.order55 import balanced_squares
OUT=ROOT/"certificates/order55_sign";exe=ROOT/"native/order55_sign/small_regression.exe"
subprocess.run(["clang++","-O3","-std=c++20","-march=native",str(ROOT/"native/order55_sign/small_regression.cpp"),"-o",str(exe)],check=True)
for order in (15,21):
 subprocess.run([str(exe),"brute",str(order),str(OUT/f"small{order}_brute.json")],check=True)
reports=[]
for n in (15,21):
 raw=json.loads((OUT/f"small{n}_brute.json").read_text());M=raw["normalized"];h=n//2
 assert 2*isqrt(n**n)+2<2**61-1
 # Full multiset universe from the sign-specific moments, independently recursive.
 records=[];weight_rows=[]
 for k in range(1,h+1):
  lower=balanced_squares(k*(k-1)//2,h);u=min(first_upper(k,n),fourth_upper(k,lower,n))
  weight_rows.append(dict(k=k,upper_floor=int(u),excluded=u<M))
  if u<M:continue
  limit=lower
  while fourth_upper(k,limit+1,n)>=M:limit+=1
  def generate(left,last,total,square,values):
   if not left:
    if total==0:records.append((k,values))
    return
   if total<last*left or square+balanced_squares(total,left)>limit:return
   for value in range(last,min(k,total//left)+1):generate(left-1,value,total-value,square+value*value,values+[value])
  generate(h,0,k*(k-1)//2,0,[])
  S2=F(n*(k*k+2*limit)-k**4,2)
  assert (2*(n-2*k))**2*S2**h < (2**61-1)**2*h**h
 path=OUT/f"small{n}_parts.txt";path.write_text("\n".join(" ".join(map(str,[k]+vs)) for k,vs in records)+"\n")
 targets=OUT/f"small{n}_profiles.txt";lifts=OUT/f"small{n}_lifts.txt"
 a=subprocess.run([str(exe),"profiles",str(n),str(path),str(targets),str(M)],text=True,capture_output=True,check=True)
 b=subprocess.run([str(exe),"lifts",str(n),str(targets),str(lifts)],text=True,capture_output=True,check=True)
 realized=[l.split() for l in lifts.read_text().splitlines()];words={w for w,d in realized}
 expected={w for w in raw["winners"] if w.count("1")<=h}
 assert words==expected
 assert max(int(d) for w,d in realized)==M
 # Every upper-half winner is exactly a complement of a lower-half winner.
 allw=words|{"".join(str(1-int(v)) for v in w) for w in words}
 assert allw==set(raw["winners"])
 assert all(raw["max_by_weight"][k]==raw["max_by_weight"][n-k] for k in range(n+1))
 # Distinct maximizing orbits are checked by integer Bareiss.
 units=[u for u in range(1,n) if __import__("math").gcd(u,n)==1]
 canonical=lambda w:min("".join(w[(u*i+t)%n] for i in range(n)) for u in units for t in range(n))
 classes={canonical(w) for w in words}
 for word in classes:
  assert abs(bareiss(circulant([2*int(v)-1 for v in word])))==raw["raw"]
 reports.append(dict(n=n,status="PASS",raw=raw["raw"],normalized=M,all_sign_words=raw["all_words"],lower_weight_winners=len(words),total_winners=len(allw),lower_affine_classes=len(classes),weights=sorted({w.count("1") for w in allw}),weight_reduction="PASS",profile_coverage="PASS",lift_coverage="PASS",winner_set_equality="PASS",formal_multisets=len(records),profile_log=a.stderr,lift_log=b.stderr,weight_screen=weight_rows))
 print(reports[-1],flush=True)
(OUT/"small_order_regression.json").write_text(json.dumps(dict(status="PASS",scope="Generic exact sign moment/profile/CRT-column-split pipeline at 15 and 21; independent direct Gray sign brute force",reports=reports),indent=2)+"\n")

