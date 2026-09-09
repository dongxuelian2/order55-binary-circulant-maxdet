"""Independent native MITM coverage control on a small exact margin fiber."""
from itertools import product
from pathlib import Path
import subprocess
from maxdet.order55 import correlations,fold
ROOT=Path(__file__).resolve().parents[1]
word='1111111'+'0'*48
c=correlations(word);r=fold(word,5);s=fold(word,11)
case=ROOT/'results/order55_mitm_control.txt'
case.write_text(' '.join(map(str,list(r)+list(s)+list(c)))+'\n')
choices=[[mask for mask in range(32) if mask.bit_count()==v] for v in s]
expected=set();assignments=0
for cols in product(*choices):
    assignments+=1
    if any(sum(mask>>i&1 for mask in cols)!=r[i] for i in range(5)):continue
    a=sum(1<<((11*i+45*j)%55) for j,mask in enumerate(cols) for i in range(5) if mask>>i&1)
    if all((a&(((a<<t)|(a>>(55-t)))&((1<<55)-1))).bit_count()==c[t] for t in range(1,28)):
        expected.add(''.join(str(a>>i&1) for i in range(55)))
assert word in expected
for split in (4,5,6):
    text=subprocess.check_output([str(ROOT/'native/order55_torus_lift.exe'),str(case),str(split),str(ROOT/f'results/order55_mitm_control_split{split}.json')],text=True)
    actual={line.split()[1] for line in text.splitlines()}
    assert actual==expected,(split,len(actual),len(expected))
print('MITM exact fiber control PASS',assignments,'assignments',len(expected),'solutions; splits 4,5,6 agree')
