"""Publication transcription check for the exact sign theorem."""
from pathlib import Path
import json
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/"certificates/order55_sign"
g=json.loads((OUT/"global.json").read_text());w=json.loads((OUT/"winner.json").read_text())
c=json.loads((OUT/"classification.json").read_text())
tex=(ROOT/"paper/order55_global/arxiv/main.tex").read_text(encoding="utf8")
needles=[g["maximum_normalized"],g["maximum_raw"],w["binary_determinant"],*w["cyclotomic_norms"],
*[r["canonical_binary"] for r in c["classes"]],w["sign_word"],"Qichao Wang","Daoyu Dong",
r"\label{thm:sign}",r"\label{lem:sign-bridge}","independent","FKM"]
needles += [f'{g[k]:,}' for k in ("formal_correlation_profiles","explicit_profile_evaluations","correlation_targets","lift_tasks","joined_words_per_full_lift")]
for needle in needles:
 if needle not in tex:raise RuntimeError("Missing sign transcription: "+needle)
for seq,key in (("A215723","maximum_raw"),("A215897","maximum_normalized")):
 draft=(ROOT/f"paper/order55_global/oeis_{seq}_draft.md").read_text()
 if g[key] not in draft:raise RuntimeError("Wrong OEIS draft")
log=(ROOT/"tmp/pdfs/joint/main.log").read_text(errors="replace")
if "Overfull" in log or "undefined references" in log:raise RuntimeError("Unresolved PDF typesetting issue")
print("ORDER-55 JOINT MANUSCRIPT AND OEIS TRANSCRIPTION: PASS")
