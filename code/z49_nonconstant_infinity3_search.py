#!/usr/bin/env python3
"""Adapted Hensel search on the 3-integral q=infinity order-49 chart.

Use the sign quotient w=u^2 and tau=v/u.  At the zero of B, write the
necessary square as t=1+tau*s, so

    r=1+2*s+tau*(s^2+2*s)+tau^2*s^2.

Two strict transforms isolate the only integral 3-adic disk.  Fixing its
parameter B leaves an invertible 4-by-4 Hensel system in the adapted
coordinates A,W,L,S below.  Reconstruction in these coordinates sees
points that may have large original coefficients.
"""

import argparse
from fractions import Fraction
from math import gcd,isqrt
import resource
import sys

import sympy as sp

sys.path.insert(0,"code")
import z49_structural_contact_iterate as contact
from z49_structural_3adic import (
    mat_vec,matrix_inverse_mod,rational_reconstruction_bounded,
)

JACOBIAN=[[1,0,0,0],[1,-1,1,0],[-1,1,1,-1],[1,1,-1,-1]]
INVERSE=matrix_inverse_mod(JACOBIAN,3)
assert INVERSE is not None


def quotient_terms():
    _x,variables,_h,_f,_solved,equations=contact.derive_system()
    a,b,u,v,r=variables;tau=sp.symbols("tau")
    output=[];maxima=[0]*5
    for equation in equations:
        poly=sp.Poly(equation.as_expr().subs(v,tau*u),a,b,u,tau,r)
        terms=[]
        for monomial,coefficient in poly.terms():
            assert monomial[2]%2==0
            qmon=(monomial[0],monomial[1],monomial[2]//2,
                  monomial[3],monomial[4])
            terms.append((qmon,int(coefficient)))
            maxima=[max(maxima[i],qmon[i]) for i in range(5)]
        output.append(terms)
    return output,maxima


def evaluate_mod(poly,values,modulus,maxima):
    powers=[]
    for value,maximum in zip(values,maxima):
        row=[1]
        for _ in range(maximum):row.append(row[-1]*value%modulus)
        powers.append(row)
    answer=0
    for monomial,coefficient in poly:
        term=coefficient
        for i,exponent in enumerate(monomial):
            if exponent:term=term*powers[i][exponent]%modulus
        answer+=term
    return answer%modulus


def evaluate_exact(poly,values):
    answer=Fraction(0)
    for monomial,coefficient in poly:
        term=Fraction(coefficient)
        for value,exponent in zip(values,monomial):term*=value**exponent
        answer+=term
    return answer


def coordinates(B,A,W,L,S,modulus=None):
    a=-8-36*B-36*B**2-27*B**3+81*A
    b=1+3*B
    w=28+18*B-27*B**2-27*B**3+81*W
    tau=18+27*B+81*L
    s=22-27*B-27*B**2+81*S
    r=1+2*s+tau*(s**2+2*s)+tau**2*s**2
    values=(a,b,w,tau,r)
    return values if modulus is None else tuple(value%modulus for value in values)


def initial(B):
    return [B**4+B**3+B**2-B+1,
            -B**2+B+1,
            B**2-B+1,
            B**2-B+1]


def strict_values(polys,maxima,B,current,modulus):
    scale=81;wide=scale*modulus
    values=coordinates(B,*current,modulus=wide)
    raw=[evaluate_mod(poly,values,wide,maxima) for poly in polys]
    assert all(value%scale==0 for value in raw)
    return [(value//scale)%modulus for value in raw]


def lift(polys,maxima,bfrac,precision):
    numerator,denominator=bfrac
    B=numerator*pow(denominator,-1,3)%3
    current=[value%3 for value in initial(B)];modulus=3
    for _ in range(1,precision):
        next_modulus=3*modulus;wide=81*next_modulus
        B=numerator*pow(denominator,-1,wide)%wide
        values=strict_values(polys,maxima,B,current,next_modulus)
        assert all(value%modulus==0 for value in values)
        rhs=[-(value//modulus)%3 for value in values]
        correction=mat_vec(INVERSE,rhs,3)
        current=[(current[i]+modulus*correction[i])%next_modulus
                 for i in range(4)]
        modulus=next_modulus
    return current,modulus


def rationals(height):
    for denominator in range(1,height+1):
        if denominator%3==0:continue
        for numerator in range(-height,height+1):
            if gcd(abs(numerator),denominator)==1:
                yield numerator,denominator


def rational_square(value):
    if value<0:return None
    numerator=isqrt(value.numerator);denominator=isqrt(value.denominator)
    if numerator*numerator!=value.numerator or denominator*denominator!=value.denominator:
        return None
    return Fraction(numerator,denominator)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--parameter-height",type=int,default=80)
    parser.add_argument("--coordinate-bound",type=int,default=5000)
    parser.add_argument("--precision",type=int,default=17)
    parser.add_argument("--memory-mb",type=int,default=300)
    args=parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,
                       (args.memory_mb*1024**2,args.memory_mb*1024**2))
    assert 3**args.precision>2*args.coordinate_bound**2
    polys,maxima=quotient_terms()
    lifted=reconstructed=exact_quotient=square_scale=exact_hits=0
    for bfrac in rationals(args.parameter_height):
        residues,modulus=lift(polys,maxima,bfrac,args.precision);lifted+=1
        rec=[rational_reconstruction_bounded(value,modulus,
                                              args.coordinate_bound)
             for value in residues]
        if any(value is None for value in rec):continue
        reconstructed+=1
        B=Fraction(*bfrac);adapted=[Fraction(*value) for value in rec]
        values=coordinates(B,*adapted)
        if any(evaluate_exact(poly,values)!=0 for poly in polys):continue
        exact_quotient+=1
        root=rational_square(values[2])
        if root is None:continue
        square_scale+=1
        a,b,w,tau,r=values
        if r==1 or a+b==Fraction(5,2) or 1+tau==0 or 1+tau*r==0:
            continue
        exact_hits+=1
        print("EXACT_OPEN_HIT","B",B,"adapted",adapted,
              "quotient",values,"u",root,"v",tau*root)
    print("Z49_NONCONSTANT_INFINITY3_SEARCH")
    print("parameter_height",args.parameter_height,
          "coordinate_bound",args.coordinate_bound,"precision",args.precision,
          "modulus",3**args.precision)
    print("lifts",lifted,"reconstructed",reconstructed,
          "exact_quotient",exact_quotient,"square_scale",square_scale,
          "exact_open_hits",exact_hits)
    print("Z49_NONCONSTANT_INFINITY3_SEARCH_DONE")


if __name__=="__main__":main()
