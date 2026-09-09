"""Exact moment product bound with a uniform spectral cap."""
from fractions import Fraction as F
from maxdet.order55 import sqrt_interval

def boxed(total,squares,count,cap):
    total,squares,cap=map(F,(total,squares,cap))
    squares=max(squares,total*total/count)
    best=F(0)
    for endpoints in range(count):
        h=count-endpoints;T=total-endpoints*cap;L=squares-endpoints*cap*cap
        if T<0 or T>h*cap or L<T*T/h or L>T*cap:continue
        mean=T/h;V=L-T*T/h
        if V==0:
            best=max(best,cap**endpoints*mean**h);continue
        for low_count in range(1,h):
            a2=V*F(h-low_count,h*low_count);b2=V*F(low_count,h*(h-low_count))
            if a2>=mean*mean:continue
            blo,bhi=sqrt_interval(b2)
            if mean+blo>cap:continue
            alo,ahi=sqrt_interval(a2)
            best=max(best,cap**endpoints*(mean-alo)**low_count*(mean+bhi)**(h-low_count))
    return best
