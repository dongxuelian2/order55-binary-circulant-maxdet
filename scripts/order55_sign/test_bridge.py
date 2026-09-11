from pathlib import Path
import sys,random,json
ROOT=Path(__file__).resolve().parents[2];sys.path.insert(0,str(ROOT/"src"))
from maxdet.sign55 import *
from maxdet.order55 import modular_determinant
count=0
for n in (3,5,7,9):
 for mask in range(1<<n):
  a=[mask>>i&1 for i in range(n)];x=[2*b-1 for b in a];k=sum(a)
  da=bareiss(circulant(a));dx=bareiss(circulant(x))
  assert k*dx==2**(n-1)*(2*k-n)*da
  assert bareiss(circulant([-v for v in x]))==-dx
  if k:assert abs(dx)*k==2**(n-1)*abs(n-2*k)*da
  assert dx%2**(n-1)==0
  c=[sum(a[i]*a[(i+s)%n] for i in range(n)) for s in range(n)]
  rho=[sum(x[i]*x[(i+s)%n] for i in range(n)) for s in range(n)]
  assert rho==[n-4*k+4*v for v in c]
  assert sum(rho)==(2*k-n)**2
  assert all(v%4==n%4 for v in rho)
  count+=1
r=random.Random(551109)
for k in range(1,28):
 a=[1]*k+[0]*(55-k);r.shuffle(a)
 da=bareiss(circulant(a));dx=bareiss(circulant([2*v-1 for v in a]))
 assert k*dx==2**54*(2*k-55)*da
 c=correlations("".join(map(str,a)))
 total=F(k*(55-k),2);square=F(55*sum(v*v for v in c)-k**4,2)
 # Validate Fourier first and fourth moments modulo a prime by independent evaluation.
 p=2305843009213696591
 z=next(pow(b,(p-1)//55,p) for b in range(2,200) if pow(pow(b,(p-1)//55,p),5,p)!=1 and pow(pow(b,(p-1)//55,p),11,p)!=1)
 q=[]
 for j in range(1,28):
  f=sum(a[t]*pow(z,j*t,p) for t in range(55))%p
  g=sum(a[t]*pow(z,(-j*t)%55,p) for t in range(55))%p
  q.append(f*g%p)
 assert (sum(q)*total.denominator-total.numerator)%p==0
 assert (sum(v*v for v in q)*square.denominator-square.numerator)%p==0
 assert abs(dx)//2**54<=fourth_upper(k,sum(v*v for v in c[1:28]))
 count+=1
old=json.loads((ROOT/"certificates/order55_fixed_weight_1_7.json").read_text())
transfer=[]
for row in old["results"]:
 k=row["weight"];d=int(row["absolute_determinant"])
 assert bareiss(circulant(list(map(int,row["word"]))))==d
 numerator=(55-2*k)*d
 assert numerator%k==0
 transfer.append(dict(k=k,normalized=str(numerator//k),source_binary_determinant=str(d)))
OUT=ROOT/"certificates/order55_sign"
(OUT/"fixed_weight_1_7.json").write_text(json.dumps(dict(status="EXACT TRANSFER OF EXISTING FIXED-WEIGHT CERTIFICATES",source="certificates/order55_fixed_weight_1_7.json",results=transfer),indent=2)+"\n")
(OUT/"bridge_regression.json").write_text(json.dumps(dict(status="PASS",exhaustive_orders=[3,5,7,9],random_order55_weights=list(range(1,28)),words_tested=count,checks=["direct sign Bareiss vs binary bridge","negation","divisibility","autocorrelation congruence","modular Fourier first/fourth moments","sign fourth bound"]),indent=2)+"\n")
print("BRIDGE AND MOMENT REGRESSION PASS",count)

