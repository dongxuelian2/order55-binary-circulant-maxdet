"""Rational Taylor enclosures and integer overflow audit for native pruning."""
from fractions import Fraction as F
from math import factorial,isqrt
from pathlib import Path
import json
S=1<<24
lo=F('3.14159265358979323846264338327950288419716939937510')
hi=lo+F(1,10**50)
rows=[]
for r in range(28):
    x=(lo+hi)*r/55
    poly=sum((-1)**j*x**(2*j)/factorial(2*j) for j in range(41))
    error=x**82/factorial(82)+(hi-lo)*r/55
    a,b=2*(poly-error)*S,2*(poly+error)*S
    rows.append((a.numerator//a.denominator,-(-b.numerator//b.denominator)))
root=Path(__file__).resolve().parents[1]
(root/'results/order55_cosine_intervals.txt').write_text('\n'.join(f'{a} {b}' for a,b in rows)+'\n')
for k in range(24,28):
    parts=json.loads((root/f'results/order55_correlation_partitions_k{k}.json').read_text())
    s2=F(55*(k*k+2*parts['maximum_half_squares'])-k**4,2)
    abs_sum=isqrt(int(27*s2))+1
    max_deviation=max(sum(abs(v-k*(k-1)//54) for v in p['values']) for p in parts['partitions'])
    error=F(27*max_deviation*max(b-a for a,b in rows),S)
    A=F(abs_sum)+error
    prefix=max([F(1)]+[(A/(16*h))**h for h in range(1,28)])
    max_state=((1<<64)+27)*prefix
    assert max_state*A*S+(1<<28) < 1<<128
print('certified cosine intervals and all native u128 intermediates PASS')
