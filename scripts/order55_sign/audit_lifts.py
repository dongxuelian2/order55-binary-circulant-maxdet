"""Independent reconstruction of every target fiber and complete orbit audit."""
from pathlib import Path
import sys,json,hashlib
from collections import defaultdict,Counter
from math import gcd
from fractions import Fraction
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/"certificates/order55_sign";CLEAN=OUT/"independent"
sys.path.insert(0,str(ROOT/"src"))
from maxdet.sign55 import certify,bareiss,circulant
def write(name,obj):(OUT/name).write_text(json.dumps(obj,indent=2)+"\n")
assert json.loads((CLEAN/"bounds_report.json").read_text())["status"]=="PASS"
assert json.loads((CLEAN/"fold_report.json").read_text())["status"]=="PASS"
targets=json.loads((OUT/"correlation_targets.json").read_text())
maps={}
for k in {t["k"] for t in targets}:
 for m in (5,11):
  lookup=defaultdict(set)
  for line in (CLEAN/f"folds_m{m}_k{k}.jsonl").read_text().splitlines():
   row=json.loads(line);a=tuple(row["profile"])
   h=tuple(sum(a[i]*a[(i+s)%m] for i in range(m)) for s in range(m));lookup[h].add(a)
  maps[k,m]={h:sorted(vs) for h,vs in lookup.items()}
expected=set();payloads={}
for idx,t in enumerate(targets):
 k=t["k"];full=tuple([k]+t["correlations"]+t["correlations"][::-1])
 rs=maps[k,5][tuple(sum(full[i::5]) for i in range(5))]
 ss=maps[k,11][tuple(sum(full[i::11]) for i in range(11))]
 assert full[11]<=15 and full[22]<=15 # Packed-key range is sufficient.
 payloads[idx]=(rs,ss,full)
 for i in range(len(rs)):
  for j in range(len(ss)):expected.add(f"{idx}_{i}_{j}")
actual=set();parts={}
for part in range(8):
 tasks=set()
 for line in (OUT/f"lift_part{part}.txt").read_text().splitlines():
  task,*fields=line.split();v=tuple(map(int,fields));idx,i,j=map(int,task.split("_"))
  rs,ss,full=payloads[idx]
  assert v==rs[i]+ss[j]+full
  assert task not in actual;actual.add(task);tasks.add(task)
 parts[part]=tasks
assert actual==expected
records={};words={}
for mode in ("production","independent"):
 rows={};allwords=set()
 for part in range(8):
  rr=[json.loads(l) for l in (OUT/f"{mode}_lift_part{part}_audit.jsonl").read_text().splitlines()]
  assert len(rr)==len(parts[part]) and {r["task"] for r in rr}==parts[part]
  for r in rr:
   assert r["split"]==(5 if mode=="production" else 6)
   rows[r["task"]]=(r["joined_words"],r["solutions"])
  ww=set((OUT/f"{mode}_lift_part{part}_words.txt").read_text().splitlines())
  counted=Counter(l.split()[0] for l in ww)
  for r in rr:assert counted[r["task"]]==r["solutions"]
  allwords.update(ww)
 for entry in allwords:
  task,w=entry.split();idx,i,j=map(int,task.split("_"));rs,ss,target=payloads[idx]
  a=tuple(map(int,w));c=tuple(sum(a[t]*a[(t+s)%55] for t in range(55)) for s in range(55))
  assert c==target
  assert tuple(sum(a[t::5]) for t in range(5))==rs[i]
  assert tuple(sum(a[t::11]) for t in range(11))==ss[j]
 records[mode]=rows;words[mode]=allwords
assert records["production"]==records["independent"]
assert words["production"]==words["independent"]
write("cleanroom_lift_report.json",dict(status="PASS",independent_full_fiber_union="PASS",target_count=len(targets),expected_tasks=len(expected),actual_tasks=len(actual),joined_words=sum(r[0] for r in records["independent"].values()),normalized_solutions=len(words["independent"]),exact_fiber_counts_equal=True,full_fiber_word_sets_equal=True,independently_reconstructed_margins=True,production_split=5,independent_split=6,proof="Every binary word can be mapped to a canonical correlation target by a unit. The independent necklace fold lists include every rotation class in that orientation. The two fold rotations lift simultaneously to one translation modulo 55. Each Cartesian margin fiber is fully enumerated by two distinct implementations."))
# Affine audit uses PUSHFORWARD SUPPORT SETS, not production's pullback words.
def orbit(w):
 support={i for i,b in enumerate(w) if b=="1"};result=set();fix=[]
 for u in range(55):
  if gcd(u,55)!=1:continue
  for t in range(55):
   transformed={(u*i+t)%55 for i in support}
   v="".join("1" if i in transformed else "0" for i in range(55))
   result.add(v)
   if transformed==support:fix.append((u,t))
 assert len(result)*len(fix)==2200
 return result,fix
certs=[certify(entry.split()[1]) for entry in words["independent"]]
maximum=max(int(c["normalized"]) for c in certs);classes=[];all_winners=set();seen=set()
for c in certs:
 if int(c["normalized"])!=maximum or c["word"] in seen:continue
 positive,fix=orbit(c["word"]);seen.update(positive)
 complement="".join(str(1-int(v)) for v in c["word"]);negative,negfix=orbit(complement)
 assert positive.isdisjoint(negative)
 for seq,orb,stabilizer in ((min(positive),positive,fix),(min(negative),negative,negfix)):
  x=[2*int(b)-1 for b in seq]
  assert abs(bareiss(circulant(x)))==maximum*2**54
  classes.append(dict(canonical_binary=seq,canonical_sign="".join("+" if v=="1" else "-" for v in seq),weight=seq.count("1"),row_sum=2*seq.count("1")-55,orbit_size=len(orb),stabilizer_size=len(stabilizer),stabilizer=stabilizer))
  all_winners.update(orb)
winner=certs[0];winner["status"]="EXACT CANDIDATE; independent profile audit required before global closure"
write("winner.json",winner)
low_binary=int(winner["binary_determinant"]);old=134694094094758395331307111329132
comparison=dict(binary_weight=winner["weight"],binary_determinant=str(low_binary),binary_complement_determinant=str(low_binary*(55-winner["weight"])//winner["weight"]),existing_binary_global=str(old),shortfall=str(old-low_binary),fraction_of_existing=str(Fraction(low_binary,old)),old_binary_induced_normalized="4810503360527085547546682547469",sign_improvement_ratio=str(Fraction(maximum,4810503360527085547546682547469)),same_as_old_binary_winner=False,half_binary_correlation_histogram=dict(Counter(winner["correlations"][1:28])),half_sign_correlation_histogram=dict(Counter(winner["sign_correlations"][1:28])))
write("comparison.json",comparison)
wordfile=OUT/"all_maximizing_sign_words.txt";wordfile.write_text("\n".join(sorted(all_winners))+"\n")
write("classification.json",dict(status="COMPLETE CLASSIFICATION OF ENUMERATED MAXIMUM; global profile audit required",maximum_normalized=str(maximum),total_maximizing_sign_words=len(all_winners),affine_orbits=len(classes),affine_plus_negation_classes=len(classes)//2,classes=classes,all_words_sha256=hashlib.sha256(wordfile.read_bytes()).hexdigest(),canonicalization="Independent pushforward support action; sign reversal creates a separate affine orbit"))
print("FULL INDEPENDENT FIBER UNION PASS",len(expected),"TASKS",flush=True)
print("CLASSIFICATION",len(all_winners),"WORDS",len(classes),"AFFINE ORBITS",len(classes)//2,"WITH NEGATION",flush=True)

