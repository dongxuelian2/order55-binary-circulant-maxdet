"""Close the sign certificate only when every required exhaustive audit passed."""
from pathlib import Path
import json,hashlib
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/"certificates/order55_sign"
if not __debug__:raise RuntimeError("Verification must run with Python assertions enabled")
def read(n):return json.loads((OUT/n).read_text())
def write(n,x):(OUT/n).write_text(json.dumps(x,indent=2)+"\n")
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
for name in ("bridge_regression.json","small_order_regression.json","cleanroom_profile_report.json","cleanroom_lift_report.json"):
 assert read(name)["status"]=="PASS",name
profiles=read("profile_coverage.json");clean=read("cleanroom_profile_report.json");lifts=read("cleanroom_lift_report.json")
assert profiles["surviving_profiles"]==clean["survivors"]==lifts["target_count"]==17835
assert lifts["expected_tasks"]==lifts["actual_tasks"]==151772
assert lifts["independent_full_fiber_union"]=="PASS"
winner=read("winner.json");classification=read("classification.json")
M=int(winner["normalized"]);raw=int(winner["raw_absolute"])
assert raw==2**54*M and M==int(classification["maximum_normalized"])
assert M==int(read("incumbent.json")["normalized"])
winner["status"]="GLOBAL MAXIMUM PROVED";classification["status"]="COMPLETE"
write("winner.json",winner);write("classification.json",classification)
pm=read("profile_manifest.json");pm["status"]="COMPLETE; production and independent universes audited";write("profile_manifest.json",pm)
lm=read("lift_manifest.json");lm["status"]="COMPLETE; every listed fiber independently reconstructed and exhaustively enumerated twice";write("lift_manifest.json",lm)
profiles["status"]="COMPLETE; exact target equality after certified exclusion of indefinite formal profiles";write("profile_coverage.json",profiles)
sources=list((ROOT/"scripts/order55_sign").glob("*.py"))+list((ROOT/"native/order55_sign").glob("*.cpp"))+[
 ROOT/"scripts/verify_order55_sign.py",ROOT/"research/order55_sign/verify_binary_companion.py",ROOT/"src/maxdet/sign55.py",ROOT/"src/maxdet/order55.py",ROOT/"src/maxdet/boxed_moment.py",
 ROOT/"scripts/build_order55_global.py",ROOT/"native/order55_profiles.cpp",ROOT/"native/order55_torus_lift.cpp",
 ROOT/"native/order55_correlation_screen.cpp",ROOT/"audit/order55/independent_torus_lift.cpp",ROOT/"audit/order55/independent_profile_screen.cpp"]
write("source_hashes.json",{p.relative_to(ROOT).as_posix():sha(p) for p in sorted(set(sources))})
assets={p.relative_to(OUT).as_posix():sha(p) for p in sorted(OUT.rglob("*")) if p.is_file() and p.name not in ("global.json","hash_manifest.json")}
write("hash_manifest.json",assets)
write("global.json",dict(status="COMPLETE",order=55,alphabet=[-1,1],maximum_normalized=str(M),maximum_raw=str(raw),
 exact_global_upper_bound_normalized=str(M),global_completeness="PASS",cleanroom_profile_audit="PASS",cleanroom_lift_audit="PASS",
 small_order_regression="PASS",winning_weights=[c["weight"] for c in classification["classes"]],winning_row_sums=[c["row_sum"] for c in classification["classes"]],
 maximizing_sign_words=classification["total_maximizing_sign_words"],affine_orbits=classification["affine_orbits"],
 affine_plus_negation_classes=classification["affine_plus_negation_classes"],screened_weights=list(range(1,28)),
 moment_active_weights=[22,23,24,25,26],weights_with_lift_targets=[23,24,25,26],
 formal_correlation_profiles=profiles["total_formal_profiles"],explicit_profile_evaluations=profiles["explicitly_enumerated"],
 boxed_pruned_profiles=profiles["boxed_pruned"],correlation_targets=clean["survivors"],independent_indefinite_exclusions=clean["indefinite_profiles_excluded"],
 lift_tasks=lifts["actual_tasks"],joined_words_per_full_lift=lifts["joined_words"],
 old_binary_global_theorem_used_for_sign_upper_bound=False,hash_manifest_sha256=sha(OUT/"hash_manifest.json")))
print("SIGN GLOBAL CERTIFICATE CLOSED",M,raw,flush=True)

