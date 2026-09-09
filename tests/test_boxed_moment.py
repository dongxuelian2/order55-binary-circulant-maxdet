from fractions import Fraction
from itertools import product
from math import prod
from maxdet.boxed_moment import boxed

def test_boxed_moment_exhaustive_grid():
    for n in (2,3,4):
        for values in product(range(1,5),repeat=n):
            for cap in (Fraction(max(values)),Fraction(max(values))+Fraction(1,10)):
                assert boxed(sum(values),sum(x*x for x in values),n,cap)>=prod(values)
