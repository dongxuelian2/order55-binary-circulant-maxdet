"""Verify the complete order-55 SIGN certificate; --replay redoes all searches."""
from pathlib import Path
import argparse,json,hashlib,sys,subprocess
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"certificates/order55_sign"
sys.path.insert(0,str(ROOT/"src"))
from maxdet.sign55 import certify,weight_screen
def require(condition,message):
 if not condition:raise RuntimeError(message)
def read(n):return json.loads((OUT/n).read_text())
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(name,*args):subprocess.run([sys.executable,str(ROOT/"scripts/order55_sign"/name),*args],cwd=ROOT,check=True)
parser=argparse.ArgumentParser();parser.add_argument("--replay",action="store_true",help="recompile and rerun complete profile, lift and independent checks")
args=parser.parse_args()
require(__debug__,"Run without Python -O")
inc=read("incumbent.json");require(certify(inc["word"])["normalized"]==inc["normalized"],"Invalid exact incumbent")
if args.replay:
 run("test_bridge.py");run("prepare.py");run("enumerate_profiles.py")
 run("run_lifts.py","production");run("run_lifts.py","independent")
 run("audit_profiles.py");run("audit_lifts.py");run("small_regression.py");run("finalize.py")
g=read("global.json");require(g["status"]=="COMPLETE","Global certificate incomplete")
assets=read("hash_manifest.json")
require(sha(OUT/"hash_manifest.json")==g["hash_manifest_sha256"],"Asset manifest changed")
current={p.relative_to(OUT).as_posix() for p in OUT.rglob("*") if p.is_file() and p.name not in ("global.json","hash_manifest.json")}
require(current==set(assets),"Unlisted or missing certificate assets")
for name,digest in assets.items():require(sha(OUT/name)==digest,f"Changed certificate {name}")
for name,digest in read("source_hashes.json").items():require(sha(ROOT/name)==digest,f"Changed proof source {name}")
w=read("winner.json");exact=certify(w["word"])
for key in ("normalized","raw_absolute","word","binary_determinant","cyclotomic_norms"):
 require(w[key]==exact[key],f"Winner mismatch: {key}")
M=int(w["normalized"]);require(int(g["maximum_normalized"])==M,"Wrong global normalization")
require(int(g["maximum_raw"])==M*2**54,"Wrong raw determinant")
require(weight_screen(M)==read("weight_screen.json"),"Weight screen mismatch")
require(read("bridge_regression.json")["status"]=="PASS","Bridge failure")
require(read("small_order_regression.json")["status"]=="PASS","Small regression failure")
cp=read("cleanroom_profile_report.json");cl=read("cleanroom_lift_report.json")
require(cp["status"]==cl["status"]=="PASS","Independent audit incomplete")
require(cl["independent_full_fiber_union"]=="PASS","Independent full fiber union absent")
# Compare every explicit set, including the independently certified
# negative-eigenvalue discharges that distinguish real Gram spectra.
parts=[p for p in read("profile_manifest.json")["partitions"] if not p["excluded"]]
target_set=set();indefinite_count=0;ordered=0
cos=[tuple(map(int,l.split())) for l in (OUT/"cosine_intervals.txt").read_text().splitlines()]
def rows(path):return [json.loads(l) for l in path.read_text().splitlines() if l]
def key(r):return r["k"],tuple(r["correlations"]),int(r["absolute_profile_product"])
for p in parts:
 stem=f'corr_k{p["k"]}_part{p["partition"]}'
 a=json.loads((OUT/f"{stem}_audit.json").read_text());b=json.loads((OUT/f"independent/{stem}_audit.json").read_text())
 require(a["ordered_profiles"]==b["ordered_profiles"]==p["ordered_count"],"Profile coverage count mismatch")
 production={key(r) for r in rows(OUT/f"{stem}_targets.jsonl")}
 clean={key(r) for r in rows(OUT/f"independent/{stem}_targets.jsonl")}
 excluded=set()
 for r in json.loads((OUT/f"independent/{stem}_indefinite.json").read_text()):
  k,c,j=r["k"],r["correlations"],r["frequency"]
  q=(k<<24)+sum(v*cos[min(j*s%55,55-j*s%55)][1] for s,v in enumerate(c,1))
  require(q==r["eigenvalue_upper_numerator"] and q<0,"Invalid indefinite-profile discharge")
  excluded.add((k,tuple(c),int(r["absolute_algebraic_product"])))
 require(clean-excluded==production and excluded<=clean,"Profile target sets disagree")
 target_set|=production;ordered+=p["ordered_count"];indefinite_count+=len(excluded)
require(len(target_set)==g["correlation_targets"]==cp["survivors"],"Target union mismatch")
require(indefinite_count==g["independent_indefinite_exclusions"],"Indefinite count mismatch")
require(ordered==g["explicit_profile_evaluations"],"Ordered count mismatch")
expected=set();joined=0;solutions=0
for part in range(8):
 taskset={l.split()[0] for l in (OUT/f"lift_part{part}.txt").read_text().splitlines()}
 require(not taskset&expected,"Duplicate lift tasks");expected|=taskset
 runs={}
 for mode in ("production","independent"):
  rr=rows(OUT/f"{mode}_lift_part{part}_audit.jsonl")
  require(len(rr)==len(taskset) and {r["task"] for r in rr}==taskset,"Missing lift fiber")
  runs[mode]={r["task"]:(r["joined_words"],r["solutions"]) for r in rr}
 require(runs["production"]==runs["independent"],"Fiber counts disagree")
 a=set((OUT/f"production_lift_part{part}_words.txt").read_text().splitlines())
 b=set((OUT/f"independent_lift_part{part}_words.txt").read_text().splitlines())
 require(a==b,"Fiber word sets disagree")
 joined+=sum(x[0] for x in runs["production"].values());solutions+=len(a)
require(len(expected)==g["lift_tasks"]==cl["actual_tasks"],"Lift task total mismatch")
require(joined==g["joined_words_per_full_lift"],"Joined word total mismatch")
require(solutions==cl["normalized_solutions"],"Solution total mismatch")
classification=read("classification.json")
require(classification["status"]=="COMPLETE","Classification incomplete")
from math import gcd
orbit_words=set()
for c in classification["classes"]:
 word=c["canonical_binary"];support={i for i,v in enumerate(word) if v=="1"};orbit=set();fix=0
 for u in range(1,55):
  if gcd(u,55)!=1:continue
  for t in range(55):
   image={(u*i+t)%55 for i in support};fix+=image==support
   orbit.add("".join("1" if i in image else "0" for i in range(55)))
 require(min(orbit)==word and len(orbit)==c["orbit_size"] and fix==c["stabilizer_size"],"Orbit error")
 require(not orbit_words&orbit,"Overlapping affine classes");orbit_words|=orbit
require(orbit_words==set((OUT/"all_maximizing_sign_words.txt").read_text().splitlines()),"Complete winner word set error")
require(len(orbit_words)==g["maximizing_sign_words"],"Maximizer count mismatch")
print("ORDER-55 SIGN CIRCULANT CERTIFICATE")
print("A215897(55) =",M)
print("A215723(55) =",M*2**54)
print("GLOBAL MAXIMUM: PASS")
print("Independent profile universe and full fiber union: PASS")
print("Maximizing sign words:",len(orbit_words))

