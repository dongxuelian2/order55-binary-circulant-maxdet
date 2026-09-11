"""Exact arithmetic for SIGN circulants; no floating point proof decisions."""
from fractions import Fraction as F
from math import prod, isqrt
from maxdet.order55 import moment_product_upper, balanced_squares, cyclotomic_norms, correlations, affine_canonical

N = 55
INITIAL_BINARY = 134694094094758395331307111329132
INITIAL_NORMALIZED = 4810503360527085547546682547469
INITIAL_RAW = 86658324567737204993560693491570615564680298496
assert INITIAL_BINARY == 28 * INITIAL_NORMALIZED
assert INITIAL_RAW == 2**54 * INITIAL_NORMALIZED

def circulant(a):
    n=len(a)
    return [[a[(j-i)%n] for j in range(n)] for i in range(n)]

def bareiss(matrix):
    a=[list(row) for row in matrix]; n=len(a); previous=1; sign=1
    for j in range(n-1):
        if not a[j][j]:
            pivot=next((i for i in range(j+1,n) if a[i][j]),None)
            if pivot is None:return 0
            a[j],a[pivot]=a[pivot],a[j];sign=-sign
        pivot=a[j][j]
        for i in range(j+1,n):
            for l in range(j+1,n):
                num=a[i][l]*pivot-a[i][j]*a[j][l]
                a[i][l],rem=divmod(num,previous)
                if rem:raise ArithmeticError("Bareiss division is not exact")
            a[i][j]=0
        previous=pivot
    return sign*a[-1][-1]

def gaussian_mod(matrix,p):
    a=[[x%p for x in row] for row in matrix];n=len(a);det=1
    for j in range(n):
        i=next((i for i in range(j,n) if a[i][j]),None)
        if i is None:return 0
        if i!=j:a[i],a[j]=a[j],a[i];det=-det
        pivot=a[j][j];det=det*pivot%p;inverse=pow(pivot,-1,p)
        for i in range(j+1,n):
            factor=a[i][j]*inverse%p
            for l in range(j+1,n):a[i][l]=(a[i][l]-factor*a[j][l])%p
    return det%p

def signed_crt(residues,primes):
    value=0;modulus=1
    for r,p in zip(residues,primes):
        value+=modulus*((r-value)*pow(modulus,-1,p)%p);modulus*=p
    return value-modulus if value>modulus//2 else value

def normalized_from_binary(word):
    a=list(map(int,word));n=len(a);k=sum(a)
    if k in (0,n):return 0
    d=bareiss(circulant(a))
    q,rem=divmod(abs(n-2*k)*abs(d),k)
    assert rem==0
    return q

def first_upper(k,n=55):
    assert n%2 and 1<=k<=n//2
    return (n-2*k)*F(k*(n-k),n-1)**((n-1)//2)

def fourth_upper(k,half_squares,n=55):
    total=F(k*(n-k),2)
    squares=F(n*(k*k+2*half_squares)-k**4,2)
    return (n-2*k)*moment_product_upper(total,squares,(n-1)//2)

def weight_screen(incumbent):
    rows=[]
    for k in range(1,28):
        minimum=balanced_squares(k*(k-1)//2,27)
        first=first_upper(k);fourth=fourth_upper(k,minimum);upper=min(first,fourth)
        row=dict(k=k,coefficient=55-2*k,first_upper_num=str(first.numerator),
                 first_upper_den=str(first.denominator),first_upper_floor=int(first),
                 fourth_upper_floor=int(fourth),upper_floor=int(upper),
                 minimum_half_squares=minimum,excluded=upper<incumbent)
        if not row["excluded"]:
            limit=minimum
            while fourth_upper(k,limit+1)>=incumbent:limit+=1
            row.update(maximum_half_squares=limit,first_excluded_upper_floor=int(fourth_upper(k,limit+1)))
        rows.append(row)
    return rows

def certify(word):
    import sympy as sp
    assert len(word)==55 and set(word)<={"0","1"}
    k=word.count("1")
    if k>27:word="".join(str(1-int(b)) for b in word);k=55-k
    word=affine_canonical(word)
    a=list(map(int,word));x=[2*v-1 for v in a];matrix=circulant(x)
    direct=bareiss(matrix)
    # These primes do not depend on Fourier roots. All three are tested.
    primes=[2305843009213693951,2305843009213693921,2305843009213693907]
    assert all(sp.isprime(p) for p in primes)
    assert prod(primes)>2*isqrt(55**55)+2
    residues=[gaussian_mod(matrix,p) for p in primes]
    crt=signed_crt(residues,primes)
    norms=cyclotomic_norms(word);fourier=2**54*(2*k-55)*prod(norms)
    binary=bareiss(circulant(a))
    assert direct==crt==fourier
    assert binary==k*prod(norms)
    assert k*direct==2**54*(2*k-55)*binary
    normalized,rem=divmod(abs(direct),2**54);assert rem==0
    c=correlations(word);rho=[sum(x[t]*x[(t+s)%55] for t in range(55)) for s in range(55)]
    assert rho==[55-4*k+4*v for v in c]
    assert 55+2*sum(rho[1:28])==(2*k-55)**2
    return dict(status="EXACT WITNESS ONLY; global maximum unproved",word=word,
                sign_word="".join("+" if b=="1" else "-" for b in word),
                encoding_decimal=str(int(word,2)),weight=k,row_sum=2*k-55,
                normalized=str(normalized),raw_absolute=str(abs(direct)),
                signed_determinant=str(direct),binary_determinant=str(binary),
                cyclotomic_norms=list(map(str,norms)),correlations=list(c),sign_correlations=rho,
                gaussian_primes=primes,gaussian_residues=residues,
                verification={"integer_Bareiss":"PASS","modular_Gaussian_CRT":"PASS","cyclotomic_resultants":"PASS"})

