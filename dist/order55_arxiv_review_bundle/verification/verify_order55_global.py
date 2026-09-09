"""Order-55 certificate audit with deterministic production-generator replay.

No heuristic search or production certificate-builder module is imported.
Use --full-lifts to rerun every lift with the same native lift executable
under the different 6+5 column split.  This is a split-dependent cross-check,
not an independently implemented full fiber-union proof.
"""
from fractions import Fraction as F
from math import isqrt,gcd,factorial,prod
from collections import Counter,defaultdict
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse,hashlib,json,subprocess,tempfile
import sympy as sp
ROOT=Path(__file__).resolve().parents[1]
P=(2305843009213696591,2305843009213697141)

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text())
def lines(p):return [json.loads(x) for x in p.read_text().splitlines() if x]
def roots(x):
    x=F(x);s=10**50;r=isqrt(x.numerator*s*s//x.denominator)
    return F(r,s),F(r if r*r*x.denominator==x.numerator*s*s else r+1,s)
def balanced(t,n):
    q,r=divmod(t,n);return (n-r)*q*q+r*(q+1)**2

def product_bound(total,squares,n,cap=None):
    total,squares=F(total),F(squares);squares=max(squares,total*total/n)
    best=F(0)
    for endpoints in range(n if cap is not None else 1):
        h=n-endpoints;T=total-(endpoints*cap if cap is not None else 0);L=squares-(endpoints*cap*cap if cap is not None else 0)
        if T<0 or L<T*T/h or L>T*T:continue
        if cap is not None and (T>h*cap or L>T*cap):continue
        mean=T/h;V=L-T*T/h;factor=cap**endpoints if cap is not None else 1
        if V==0:best=max(best,factor*mean**h);continue
        for m in range(1,h):
            va=V*F(h-m,h*m);vb=V*F(m,h*(h-m))
            if va>=mean*mean:continue
            al,ah=roots(va);bl,bh=roots(vb)
            if cap is not None and mean+bl>cap:continue
            best=max(best,factor*(mean-al)**m*(mean+bh)**(h-m))
    return best

def fourth(k,squares):return (55-k)*product_bound(F(k*(55-k),2),F(55*(k*k+2*squares)-k**4,2),27)
def canonical(w):return min(''.join(w[(u*i+t)%55] for i in range(55)) for u in range(1,55) if gcd(u,55)==1 for t in range(55))
def corr(w):return [sum(int(w[i])*int(w[(i+s)%55]) for i in range(55)) for s in range(55)]
def fold(w,m):return [sum(map(int,w[i::m])) for i in range(m)]
def det(w):return int(sp.Matrix(55,55,lambda i,j:int(w[(j-i)%55])).det(method='domain-ge'))

def audit(folder,full):
    if not __debug__:raise RuntimeError('verification requires assertions enabled')
    g=load(folder/'global.json');assert g['status']=='COMPLETE' and g['order']==55
    manifest=load(folder/'hash_manifest.json');assert digest(folder/'hash_manifest.json')==g['hash_manifest_sha256'] and manifest==g['components']
    for name,h in manifest.items():assert Path(name).name==name and digest(folder/name)==h
    for name,h in load(folder/'source_hashes.json').items():assert digest(ROOT/name)==h
    M=int(g['screening_incumbent']);maximum=int(g['maximum']);assert maximum>=M
    assert all(sp.isprime(p) and (p-1)%55==0 for p in P)
    weights=load(folder/'weight_screen.json');assert [r['k'] for r in weights]==list(range(1,28))
    limits={}
    for row in weights:
        k=row['k'];b=balanced(k*(k-1)//2,27)
        U=int(min((55-k)*F(k*(55-k),54)**27,fourth(k,b)));assert row['upper_floor']==U
        assert row['excluded']==(U<M) and row['complement_weight']==55-k
        if not row['excluded']:
            t=row['maximum_half_squares'];assert int(fourth(k,t))>=M and int(fourth(k,t+1))<M;limits[k]=t
    assert sorted(limits)==g['active_weights']
    def arctan(q):
        s=sum(F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(100))
        return s,s+F(1,201*q**201)
    a,b=arctan(5);c,d=arctan(239);pilo,pihi=16*a-4*d,16*b-4*c
    cos=[tuple(map(int,l.split())) for l in (folder/'cosine_intervals.txt').read_text().splitlines()];assert len(cos)==28
    for r,(low,high) in enumerate(cos):
        x=(pilo+pihi)*r/55;y=sum((-1)**j*x**(2*j)/factorial(2*j) for j in range(43))
        error=x**86/factorial(86)+(pihi-pilo)*r/55
        assert F(low,1<<24)<=2*(y-error)<=2*(y+error)<=F(high,1<<24)
    expected=[]
    for m in (5,11):
        h=(m-1)//2;rest=27-h
        assert F(m*m*(55//m)**2,m-1)**h<P[0]
        for k in sorted(limits):
            for sq in range(k*k//m,(55//m)*k+1):
                T=F(m*sq-k*k,2);R=F(k*(55-k),2)-T
                if T<=0 or R<=0 or (55-k)*(T/h)**h*(R/rest)**rest<M:continue
                need=F(M,55-k)/(R/rest)**rest;expected.append((m,k,sq,-(-need.numerator//need.denominator)))
    assert expected==[tuple(map(int,x.split())) for x in (folder/'profile_thresholds.txt').read_text().splitlines()]
    with tempfile.TemporaryDirectory(prefix='order55-audit-') as tmp:
        tmp=Path(tmp);executables={}
        for record in load(folder/'build.json'):
            source=next(ROOT/'native'/Path(x).name for x in record['command'] if x.endswith('.cpp'))
            assert digest(source)==record['source_sha256']
            exe=tmp/(source.stem+'.exe');subprocess.run(['clang++','-O3','-std=c++20','-march=native',str(source),'-o',str(exe)],check=True);executables[source.stem]=exe
        lookup={};profile_counts=[]
        for k in sorted(limits):
            for m in (5,11):
                stored=folder/f'profiles_m{m}_k{k}.jsonl';replay=tmp/stored.name
                with replay.open('w') as f:subprocess.run([str(executables['order55_profiles']),str(folder/'profile_thresholds.txt'),str(m),str(k)],stdout=f,stderr=subprocess.DEVNULL,check=True)
                assert digest(replay)==digest(stored)
                mapping=defaultdict(set)
                for row in lines(stored):
                    a=row['profile'];assert len(a)==m and sum(a)==k and all(0<=v<=55//m for v in a)
                    x=sp.Symbol('x');N=int(sp.resultant(sum(v*x**i for i,v in enumerate(a)),sp.cyclotomic_poly(m,x),x));assert N==row['norm']
                    for u in range(1,m):
                        b=tuple(a[u*i%m] for i in range(m));b=min(b[t:]+b[:t] for t in range(m))
                        h=tuple(sum(b[i]*b[(i+s)%m] for i in range(m)) for s in range(m));mapping[h].add(b)
                lookup[k,m]=mapping
                expected_sig=sorted({h[:(m+1)//2] for h in mapping})
                assert expected_sig==[tuple(map(int,x.split())) for x in (folder/f'foldcorr_m{m}_k{k}.txt').read_text().splitlines()]
        print('HASHES, WEIGHTS, COSINE INTERVALS, FOLDED PROFILES: PASS',flush=True)
        entries=load(folder/'correlation_profiles.json');alltargets=[]
        for k in sorted(limits):
            multisets=set()
            def visit(vs,start,total,sq):
                n=27-len(vs)
                if not n:
                    if total==0:multisets.add(tuple(vs))
                    return
                if total<start*n or sq+balanced(total,n)>limits[k]:return
                for v in range(start,min(k,total//n)+1):visit(vs+[v],v,total-v,sq+v*v)
            visit([],0,k*(k-1)//2,0)
            assert multisets=={tuple(r['values']) for r in entries if r['k']==k}
        for entry in entries:
            k=entry['k'];vs=entry['values'];sq=sum(v*v for v in vs);assert sq==entry['squares']
            count=factorial(27)//prod(factorial(v) for v in Counter(vs).values());assert count==entry['ordered_count']
            base=k*(k-1)//54;ds=sorted(v-base for v in vs);caps=[]
            for j in (1,5,11):
                cs=sorted(cos[min(j*s%55,55-j*s%55)] for s in range(1,28));caps.append(F(k-base)+sum(F(d*(v[1] if d>=0 else v[0]),1<<24) for d,v in zip(ds,cs)))
            cap=max(caps);assert cap==F(entry['cap'])
            s2=F(55*(k*k+2*sq)-k**4,2);bound=(55-k)*product_bound(F(k*(55-k),2),s2,27,cap);assert int(bound)==entry['upper_floor']
            assert (2*(55-k))**2*s2**27<(P[0]*P[1])**2*27**27
            A=F(isqrt(int(27*s2))+1)+F(27*sum(abs(v) for v in ds)*max(b-a for a,b in cos),1<<24);pref=max([F(1)]+[(A/(16*h))**h for h in range(1,28)])
            assert ((1<<64)+27)*pref*A*(1<<24)+(1<<28)<1<<128
            if entry['route']=='boxed moment':assert int(bound)<M;continue
            idx=entry['partition'];source=folder/f'corr_k{k}_part{idx}.txt';stored=folder/f'corr_k{k}_part{idx}_targets.jsonl';replay=tmp/stored.name;ra=tmp/(source.stem+'_audit.json')
            assert list(map(int,source.read_text().split()))==vs
            with replay.open('w') as f:subprocess.run([str(executables['order55_correlation_screen']),str(k),str(M),str(source),str(folder/'foldcorr_m'),str(ra),str(folder/'cosine_intervals.txt')],stdout=f,check=True)
            assert digest(replay)==digest(stored);fresh=load(ra);old=entry['audit'];assert fresh['ordered_profiles']==count
            for key in ('ordered_profiles','folded_matching_profiles','canonical_determinants','above_screen','max_absolute_profile_product'):assert fresh[key]==old[key]
            for t in lines(stored):t['source_partition']=idx;alltargets.append(t)
        alltargets.sort(key=lambda t:(-int(t['absolute_profile_product']),t['k'],t['correlations']))
        coverage=load(folder/'lift_targets.json');assert [r['target'] for r in coverage]==alltargets and len(coverage)==g['correlation_targets']
        tasks={}
        for idx,row in enumerate(coverage):
            assert row['id']==idx;t=row['target'];k=t['k'];c=[k]+t['correlations']+list(reversed(t['correlations']))
            rs=sorted(lookup[k,5][tuple(sum(c[i::5]) for i in range(5))]);ss=sorted(lookup[k,11][tuple(sum(c[i::11]) for i in range(11))])
            assert [list(r) for r in rs]==row['row_profiles'] and [list(s) for s in ss]==row['column_profiles'];assert len(rs)*len(ss)==row['tasks']
            for i,r in enumerate(rs):
                for j,s in enumerate(ss):tasks[f'{idx}_{i}_{j}']=list(r)+list(s)+c
        seen=set();joined=0;wordrecords=[];parts=load(folder/'lift_coverage.json')
        for part in parts:
            i=part['part'];inp=folder/f'lift_part{i}.txt';auditpath=folder/f'lift_part{i}_audit.jsonl';words=folder/f'lift_part{i}_words.txt';expected_ids=set()
            for line in inp.read_text().splitlines():
                a=line.split();task=a[0];assert task not in seen and list(map(int,a[1:]))==tasks[task];seen.add(task);expected_ids.add(task)
            audits=lines(auditpath);assert len(audits)==len(expected_ids) and {a['task'] for a in audits}==expected_ids;assert sum(a['joined_words'] for a in audits)==part['joined_words'];joined+=part['joined_words']
            actual=words.read_text().splitlines();assert len(actual)==part['solutions']==sum(a['solutions'] for a in audits)
            for line in actual:
                task,w=line.split();assert task in expected_ids and len(w)==55 and set(w)<={'0','1'};values=tasks[task];assert fold(w,5)==values[:5] and fold(w,11)==values[5:16] and corr(w)==values[16:]
                high=''.join('1' if b=='0' else '0' for b in w);value=det(high);low=det(w);assert low*(55-w.count('1'))==value*w.count('1');assert value==int(coverage[int(task.split('_')[0])]['target']['absolute_profile_product']);wordrecords.append((high,value))
        assert set(tasks)==seen and len(tasks)==g['lift_tasks'] and joined==g['joined_words']
        print('ALL CORRELATION PROFILES AND LIFT TASK COVERAGE: PASS',flush=True)
        if full:
            def replay_lift(part):
                i=part['part'];freshwords=tmp/f'lift{i}.txt';freshaudit=tmp/f'lift{i}.jsonl'
                with freshwords.open('w') as f:subprocess.run([str(executables['order55_torus_lift']),str(folder/f'lift_part{i}.txt'),'6',str(freshaudit),'batch'],stdout=f,check=True)
                assert set(freshwords.read_text().splitlines())==set((folder/f'lift_part{i}_words.txt').read_text().splitlines());old={r['task']:r for r in lines(folder/f'lift_part{i}_audit.jsonl')};new=lines(freshaudit);assert len(new)==len(old)
                for r in new:assert (r['joined_words'],r['solutions'])==(old[r['task']]['joined_words'],old[r['task']]['solutions'])
                print('SPLIT-DEPENDENT 6+5 LIFT REPLAY PASS',i,flush=True)
            with ThreadPoolExecutor(max_workers=8) as pool:list(pool.map(replay_lift,parts))
    assert max(v for w,v in wordrecords)==maximum
    classes={canonical(w) for w,v in wordrecords if v==maximum};winner=load(folder/'winner.json');assert classes=={r['word'] for r in winner['classes']} and len(classes)==g['affine_class_count']
    count=0
    for row in winner['classes']:
        w=row['word'];assert det(w)==maximum and w.count('1')==row['weight'];assert row['residues']==[maximum%p for p in P]
        complement=''.join('1' if b=='0' else '0' for b in w);assert row['complement']==complement and det(complement)*row['weight']==maximum*(55-row['weight'])
        x=sp.Symbol('x');f=sum(int(b)*x**i for i,b in enumerate(w));norms=[int(sp.resultant(f,sp.cyclotomic_poly(m,x),x)) for m in (5,11,55)];assert list(map(str,norms))==row['norms'] and w.count('1')*prod(norms)==maximum
        orbit={''.join(w[(u*i+t)%55] for i in range(55)) for u in range(1,55) if gcd(u,55)==1 for t in range(55)};assert len(orbit)==row['orbit_size'] and row['stabilizer_size']*len(orbit)==2200
        assert corr(w)==row['correlations'] and fold(w,5)==row['mod5'] and fold(w,11)==row['mod11'];count+=len(orbit)
    assert count==g['maximizing_words']==winner['total_maximizing_words']
    print('ORDER 55 GLOBAL CERTIFICATE');print('winner =',winner['classes'][0]['word']);print('weight =',winner['classes'][0]['weight']);print('determinant =',maximum)
    print('affine_class_count =',len(classes));print('maximizing_words =',count);print('domain_size =',1 << 55,'(mathematical domain cardinality; not a derived coverage count)');print('GLOBAL MAXIMUM: PASS')

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('certificate',nargs='?',type=Path,default=ROOT/'certificates/order55_global');p.add_argument('--full-lifts',action='store_true');a=p.parse_args();audit(a.certificate.resolve(),a.full_lifts)
