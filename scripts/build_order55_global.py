"""Reproduce the order-55 global certificate from exact arithmetic and C++20.

The initial witness fixes a proof threshold; no search score is trusted.
"""
from concurrent.futures import ThreadPoolExecutor, as_completed
from fractions import Fraction as F
from collections import Counter,defaultdict
from math import factorial,isqrt,gcd
from pathlib import Path
import argparse,hashlib,json,subprocess,time
import sympy as sp
from maxdet.order55 import (PRIMES,exact55,cyclotomic_norms,affine_canonical,
    correlations,fold,sector_moments,balanced_squares,fourth_bound,ryser_global)
from maxdet.boxed_moment import boxed
ROOT=Path(__file__).resolve().parents[1]
WITNESS='1111010001100011100000101010001101110110111010010010110'
M=134694094094758395331307111329132

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(path,data):path.write_text(json.dumps(data,indent=2)+'\n',encoding='utf8')
def readlines(path):return [json.loads(x) for x in path.read_text(encoding='utf-8-sig').splitlines() if x]

def cosine_table():
    # Machin's identity gives pi with certified alternating arctangent series.
    def atan_interval(q,n=90):
        a=sum(F((-1)**j,(2*j+1)*q**(2*j+1)) for j in range(n))
        b=a+F((-1)**n,(2*n+1)*q**(2*n+1))
        return min(a,b),max(a,b)
    a,b=atan_interval(5);c,d=atan_interval(239)
    lo,hi=16*a-4*d,16*b-4*c
    result=[];S=1<<24
    for r in range(28):
        x=(lo+hi)*r/55
        poly=sum((-1)**j*x**(2*j)/factorial(2*j) for j in range(41))
        error=x**82/factorial(82)+(hi-lo)*r/55
        low,high=2*(poly-error)*S,2*(poly+error)*S
        result.append((low.numerator//low.denominator,-(-high.numerator//high.denominator)))
    return result

def partitions(k,limit):
    out=[]
    def rec(lo,n,total,squares,vs):
        if not n:
            if total==0:
                count=factorial(27)
                for v in Counter(vs).values():count//=factorial(v)
                out.append(dict(values=vs,squares=squares,ordered_count=count))
            return
        if total<lo*n or squares+balanced_squares(total,n)>limit:return
        for v in range(lo,min(k,total//n)+1):rec(v,n-1,total-v,squares+v*v,vs+[v])
    rec(0,27,k*(k-1)//2,0,[])
    return out

def bound_partition(k,p,cos):
    base=k*(k-1)//54;ds=sorted(v-base for v in p['values']);caps=[]
    for j in (1,5,11):
        cs=sorted(cos[min(j*s%55,55-j*s%55)] for s in range(1,28))
        caps.append(F(k-base)+sum(F(v*(pair[1] if v>=0 else pair[0]),1<<24) for v,pair in zip(ds,cs)))
    cap=max(caps);s1=F(k*(55-k),2);s2=F(55*(k*k+2*p['squares'])-k**4,2)
    # Independent signed CRT capacity and dyadic-product overflow inequalities.
    assert (2*(55-k))**2*s2**27 < (PRIMES[0]*PRIMES[1])**2*27**27
    err=F(27*sum(abs(v) for v in ds)*max(b-a for a,b in cos),1<<24)
    A=F(isqrt(int(27*s2))+1)+err
    prefix=max([F(1)]+[(A/(16*h))**h for h in range(1,28)])
    assert ((1<<64)+27)*prefix*A*(1<<24)+(1<<28)<1<<128
    return dict(cap=str(cap),upper_floor=int((55-k)*boxed(s1,s2,27,cap)))

def run_native(name,args,output):
    cmd=[str(ROOT/f'native/{name}.exe')]+list(map(str,args))
    start=time.time()
    with output.open('w') as f:
        result=subprocess.run(cmd,cwd=ROOT,stdout=f,stderr=subprocess.PIPE,text=True,check=True)
    return dict(command=cmd,seconds=time.time()-start,stderr=result.stderr,
                source_sha256=sha(ROOT/f'native/{name}.cpp'),output_sha256=sha(output))

def build(out,workers):
    out.mkdir(parents=True,exist_ok=True)
    write(out/'global.json',dict(status='INCOMPLETE',screening_incumbent=str(M)))
    assert exact55(WITNESS)[0]==M
    compilation=[]
    for name in ('order55_profiles','order55_correlation_screen','order55_torus_lift'):
        cmd=['clang++','-O3','-std=c++20','-march=native','-Wall','-Wextra',str(ROOT/f'native/{name}.cpp'),'-o',str(ROOT/f'native/{name}.exe')]
        subprocess.run(cmd,check=True);compilation.append(dict(command=cmd,source_sha256=sha(ROOT/f'native/{name}.cpp')))
    write(out/'build.json',compilation)
    cos=cosine_table();cospath=out/'cosine_intervals.txt'
    cospath.write_text('\n'.join(f'{a} {b}' for a,b in cos)+'\n')
    weights=[];active=[];allparts={}
    for k in range(1,28):
        b=balanced_squares(k*(k-1)//2,27);u=int(min(ryser_global(k),fourth_bound(k,b)))
        row=dict(k=k,complement_weight=55-k,upper_floor=u,excluded=u<M)
        if u>=M:
            t=b
            while int(fourth_bound(k,t+1))>=M:t+=1
            row.update(minimum_half_squares=b,maximum_half_squares=t,first_excluded_upper=int(fourth_bound(k,t+1)))
            active.append(k);allparts[k]=partitions(k,t)
        weights.append(row)
    write(out/'weight_screen.json',weights)
    thresholds=[]
    for m in (5,11):
        h=(m-1)//2;rest=27-h;cap=55//m
        assert F(m*m*cap*cap,m-1)**h<PRIMES[0]
        for k in active:
            for squares in range(k*k//m,cap*k+1):
                T=F(m*squares-k*k,2);R=F(k*(55-k),2)-T
                if T<=0 or R<=0 or (55-k)*(T/h)**h*(R/rest)**rest<M:continue
                need=F(M,55-k)/(R/rest)**rest
                thresholds.append((m,k,squares,-(-need.numerator//need.denominator)))
    thresholdpath=out/'profile_thresholds.txt'
    thresholdpath.write_text('\n'.join(' '.join(map(str,row)) for row in thresholds)+'\n')
    profileaudit=[];lookups={}
    for k in active:
        for m in (5,11):
            path=out/f'profiles_m{m}_k{k}.jsonl'
            audit=run_native('order55_profiles',[thresholdpath,m,k],path)
            profiles=readlines(path);lookup=defaultdict(set)
            for row in profiles:
                a=row['profile']
                for u in range(1,m):
                    b=tuple(a[u*i%m] for i in range(m));b=min(b[t:]+b[:t] for t in range(m))
                    h=tuple(sum(b[t]*b[(t+s)%m] for t in range(m)) for s in range(m))
                    lookup[h].add(b)
            lookups[k,m]=lookup
            signatures=sorted({h[:(m+1)//2] for h in lookup})
            (out/f'foldcorr_m{m}_k{k}.txt').write_text('\n'.join(' '.join(map(str,s)) for s in signatures)+'\n')
            audit.update(k=k,m=m,canonical_survivors=len(profiles),folded_correlation_signatures=len(signatures));profileaudit.append(audit)
    write(out/'cyclotomic_profiles.json',profileaudit)
    corraudit=[];targets=[]
    for k in active:
        parts=allparts[k]
        for i,p in enumerate(parts):
            pb=bound_partition(k,p,cos);entry=dict(k=k,partition=i,**p,**pb)
            if pb['upper_floor']<M:entry['route']='boxed moment'
            else:
                entry['route']='exhaustive correlation enumeration'
                path=out/f'corr_k{k}_part{i}.txt';path.write_text(' '.join(map(str,p['values']))+'\n')
                targetpath=out/f'corr_k{k}_part{i}_targets.jsonl';auditpath=out/f'corr_k{k}_part{i}_audit.json'
                execution=run_native('order55_correlation_screen',[k,M,path,str(out/'foldcorr_m'),auditpath,cospath],targetpath)
                a=json.loads(auditpath.read_text());assert a['ordered_profiles']==p['ordered_count']
                entry.update(execution=execution,audit=a)
                for t in readlines(targetpath):
                    t['source_partition']=i;targets.append(t)
            corraudit.append(entry)
        print('CORRELATIONS COMPLETE',k,'targets',sum(t['k']==k for t in targets),flush=True)
    write(out/'correlation_profiles.json',corraudit)
    targets.sort(key=lambda t:(-int(t['absolute_profile_product']),t['k'],t['correlations']))
    taskfiles=[out/f'lift_part{i}.txt' for i in range(workers)]
    handles=[p.open('w') for p in taskfiles];coverage=[];total=0
    for idx,t in enumerate(targets):
        k=t['k'];c=[k]+t['correlations']+list(reversed(t['correlations']))
        rs=sorted(lookups[k,5][tuple(sum(c[i::5]) for i in range(5))]);ss=sorted(lookups[k,11][tuple(sum(c[i::11]) for i in range(11))])
        for i,r in enumerate(rs):
            for j,s in enumerate(ss):
                task=f'{idx}_{i}_{j}'
                handles[total%workers].write(task+' '+' '.join(map(str,list(r)+list(s)+c))+'\n');total+=1
        coverage.append(dict(id=idx,target=t,row_profiles=rs,column_profiles=ss,tasks=len(rs)*len(ss)))
    for h in handles:h.close()
    write(out/'lift_targets.json',coverage)
    print('LIFT TASKS',total,flush=True)
    def lift(i):
        auditpath=out/f'lift_part{i}_audit.jsonl';wordpath=out/f'lift_part{i}_words.txt'
        execution=run_native('order55_torus_lift',[taskfiles[i],5,auditpath,'batch'],wordpath)
        expected={line.split()[0] for line in taskfiles[i].read_text().splitlines()}
        audit=readlines(auditpath);assert len(audit)==len(expected) and {a['task'] for a in audit}==expected
        print('LIFT COMPLETE',i,len(audit),'solutions',sum(a['solutions'] for a in audit),flush=True)
        return dict(part=i,execution=execution,tasks=len(audit),joined_words=sum(a['joined_words'] for a in audit),solutions=sum(a['solutions'] for a in audit))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        liftaudit=list(pool.map(lift,range(workers)))
    assert sum(x['tasks'] for x in liftaudit)==total
    write(out/'lift_coverage.json',liftaudit)
    words=[]
    for i in range(workers):
        for line in (out/f'lift_part{i}_words.txt').read_text().splitlines():
            task,w=line.split();k=w.count('1');det,res=exact55(w)
            actual=correlations(w);target=coverage[int(task.split('_')[0])]['target']
            assert list(actual[1:28])==target['correlations']
            high=''.join('1' if b=='0' else '0' for b in w)
            highdet,highres=exact55(high);assert highdet==int(target['absolute_profile_product'])
            words.append(dict(task=task,word=high,weight=55-k,determinant=str(highdet),residues=highres))
    assert words,'the screening witness must lift'
    maximum=max(int(x['determinant']) for x in words);assert maximum>=M
    classes=sorted({affine_canonical(x['word']) for x in words if int(x['determinant'])==maximum})
    classified=[]
    for w in classes:
        orbit={''.join(w[(u*i+t)%55] for i in range(55)) for u in range(1,55) if gcd(u,55)==1 for t in range(55)}
        det,res=exact55(w);a=list(map(int,w));independent=int(sp.Matrix(55,55,lambda i,j:a[(j-i)%55]).det(method='domain-ge'));norms=cyclotomic_norms(w)
        assert det==independent==w.count('1')*norms[0]*norms[1]*norms[2]==maximum
        classified.append(dict(word=w,weight=w.count('1'),complement=''.join('1' if b=='0' else '0' for b in w),determinant=str(det),residues=res,norms=list(map(str,norms)),stabilizer_size=2200//len(orbit),orbit_size=len(orbit),correlations=correlations(w),mod5=fold(w,5),mod11=fold(w,11),sector_moments={m:list(map(str,p)) for m,p in sector_moments(w).items()}))
    write(out/'winner.json',dict(maximum=str(maximum),classes=classified,total_maximizing_words=sum(x['orbit_size'] for x in classified),verification=['finite-field Fourier CRT','SymPy domain-ge','cyclotomic resultants']))
    write(out/'lifted_words.json',words)
    source_paths=['scripts/build_order55_global.py','scripts/verify_order55_global.py','src/maxdet/order55.py','src/maxdet/boxed_moment.py']+[f'native/{name}.cpp' for name in ('order55_profiles','order55_correlation_screen','order55_torus_lift')]
    write(out/'source_hashes.json',{name:sha(ROOT/name) for name in source_paths})
    manifest={p.name:sha(p) for p in sorted(out.iterdir()) if p.is_file() and p.name not in ('global.json','hash_manifest.json')}
    write(out/'hash_manifest.json',manifest)
    write(out/'global.json',dict(status='COMPLETE',order=55,screening_incumbent=str(M),maximum=str(maximum),affine_class_count=len(classes),maximizing_words=sum(x['orbit_size'] for x in classified),active_weights=active,correlation_targets=len(targets),lift_tasks=total,joined_words=sum(x['joined_words'] for x in liftaudit),hash_manifest_sha256=sha(out/'hash_manifest.json'),components=manifest))
    print('GLOBAL CERTIFICATE GENERATED',maximum,'classes',len(classes),'words',sum(x['orbit_size'] for x in classified),flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=ROOT/'certificates/order55_global');p.add_argument('--workers',type=int,default=8);a=p.parse_args();build(a.output.resolve(),a.workers)
