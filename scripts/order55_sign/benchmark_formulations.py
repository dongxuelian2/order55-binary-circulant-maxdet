"""Representation benchmark only; timings never enter a pruning decision."""
from pathlib import Path
import json,time
OUT=Path(__file__).resolve().parents[2]/"certificates/order55_sign"
P=2305843009213696591
z=next(pow(b,(P-1)//55,P) for b in range(2,1000)
       if pow(pow(b,(P-1)//55,P),5,P)!=1 and pow(pow(b,(P-1)//55,P),11,P)!=1)
table=[[(pow(z,j*s,P)+pow(z,(-j*s)%55,P))%P for s in range(1,28)] for j in range(1,28)]
rows=json.loads((OUT/"correlation_targets.json").read_text())
sample=rows[::max(1,len(rows)//1000)]
def binary(row):
 k=row["k"];b=k*(k-1)//54;ds=[v-b for v in row["correlations"]];value=55-2*k
 for coeff in table:value=value*(k-b+sum(d*c for d,c in zip(ds,coeff) if d))%P
 return value
def sign(row):
 k=row["k"];b=k*(k-1)//54;rho=[55-4*k+4*c for c in row["correlations"]]
 baseline=55-4*k+4*b;ds=[v-baseline for v in rho];value=55-2*k
 for coeff in table:value=value*(55-baseline+sum(d*c for d,c in zip(ds,coeff) if d))%P
 return value*pow(pow(4,27,P),-1,P)%P
start=time.perf_counter();a=[binary(r) for r in sample];ta=time.perf_counter()-start
start=time.perf_counter();b=[sign(r) for r in sample];tb=time.perf_counter()-start
assert a==b
result=dict(status="PASS",scope="Representation benchmark only; not a proof dependency",profiles=len(sample),binary_seconds=ta,direct_sign_seconds=tb,exact_modular_products_equal=True,decision="Keep sparse binary c deviations: direct sign deviations are exactly four times these, with no extra independent congruence restriction at fixed weight.")
(OUT/"formulation_benchmark.json").write_text(json.dumps(result,indent=2)+"\n")
print(result)

