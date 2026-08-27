#!/usr/bin/env python3
"""Bounded fixed-slope search on the explicit order-49 incidence curve.

Fix T=v/u (u!=0), so v=T*u.  For every small rational T, use each matching
open F_5 and F_11 incidence point whose 4x4 Jacobian in (a,b,u,r) is
invertible.  Hensel lift, reconstruct those four coordinates, then exact-test
the original equations.  This avoids symbolic elimination for the much
harder nonconstant-B slices.
"""

from fractions import Fraction
from math import gcd
import argparse
import sys

sys.path.insert(0,"code")
import z49_structural_contact_iterate as contact
from z49_structural_3adic import (
    evaluate_all_mod, evaluate_exact, mat_vec, matrix_inverse_mod,
    rational_reconstruction_bounded, term_data,
)


BRANCHES = {
    5: [
        (2,1,1,3,0), (2,1,2,2,2),
        (2,1,3,3,2), (2,1,4,2,0),
    ],
    11: [
        (5,7,3,0,7), (5,7,8,0,7),
        (6,7,4,5,6), (6,7,7,6,6),
        (7,6,2,2,9), (7,6,5,4,8),
        (7,6,6,7,8), (7,6,9,9,9),
        (10,1,3,7,3), (10,1,8,4,3),
    ],
}


def rationals(height):
    out=[]
    for denominator in range(1,height+1):
        for numerator in range(-height,height+1):
            if gcd(abs(numerator),denominator)!=1:
                continue
            value=(numerator,denominator)
            if value not in out:
                out.append(value)
    return out


def lift(polys,maxima,inverse,initial,tfrac,p,precision):
    # current coordinates are a,b,u,r; v=T*u.
    current=[initial[0]%p,initial[1]%p,initial[2]%p,initial[4]%p]
    modulus=p
    numerator,denominator=tfrac
    for _ in range(1,precision):
        next_modulus=modulus*p
        tt=numerator*pow(denominator,-1,next_modulus)%next_modulus
        values=[current[0],current[1],current[2],
                tt*current[2]%next_modulus,current[3]]
        fvalues=evaluate_all_mod(polys,values,next_modulus,maxima)
        assert all(value%modulus==0 for value in fvalues)
        rhs=[-(value//modulus)%p for value in fvalues]
        correction=mat_vec(inverse,rhs,p)
        current=[(current[i]+modulus*correction[i])%next_modulus
                 for i in range(4)]
        modulus=next_modulus
    return current,modulus


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--slope-height",type=int,default=50)
    parser.add_argument("--coordinate-bound",type=int,default=5000)
    args=parser.parse_args()
    precisions={5:12,11:8}
    assert all(p**precisions[p]>2*args.coordinate_bound**2 for p in precisions)

    _x,variables,_h,_f,_solved,equations=contact.derive_system()
    polys=[term_data(poly) for poly in equations]
    maxima=[max(monomial[j] for terms in polys for monomial,_ in terms)
            for j in range(5)]
    slopes=[value for value in rationals(args.slope_height)
            if Fraction(*value)!=-1]
    attempted=singular=reconstructed=exact_hits=0
    hits=[]
    for p,branches in BRANCHES.items():
        for tfrac in slopes:
            numerator,denominator=tfrac
            if denominator%p==0:
                continue
            tres=numerator*pow(denominator,-1,p)%p
            for branch in branches:
                if branch[2]==0 or branch[3]*pow(branch[2],-1,p)%p!=tres:
                    continue
                seed={variable:value for variable,value in zip(variables,branch)}
                jacobian=[]
                for poly in equations:
                    jacobian.append([
                        int(poly.diff(variables[0]).subs(seed))%p,
                        int(poly.diff(variables[1]).subs(seed))%p,
                        int((poly.diff(variables[2])+tres*poly.diff(variables[3])).subs(seed))%p,
                        int(poly.diff(variables[4]).subs(seed))%p,
                    ])
                inverse=matrix_inverse_mod(jacobian,p)
                if inverse is None:
                    singular+=1
                    continue
                attempted+=1
                residues,modulus=lift(polys,maxima,inverse,branch,tfrac,p,
                                      precisions[p])
                rec=[rational_reconstruction_bounded(value,modulus,
                                                      args.coordinate_bound)
                     for value in residues]
                if any(value is None for value in rec):
                    continue
                reconstructed+=1
                aq,bq,uq,rq=[Fraction(*value) for value in rec]
                tq=Fraction(*tfrac);vq=tq*uq
                point=[aq,bq,uq,vq,rq]
                if not all(evaluate_exact(poly,point)==0 for poly in polys):
                    continue
                if (rq==1 or aq+bq==Fraction(5,2) or
                    uq+vq==0 or uq+vq*rq==0):
                    continue
                exact_hits+=1;hits.append(point)
                print("EXACT_OPEN_HIT","p",p,"T",tq,"point",point)

    print("Z49_STRUCTURAL_SLOPE_HENSEL")
    print("slope_height",args.slope_height,"slopes",len(slopes),
          "coordinate_bound",args.coordinate_bound)
    print("attempted_branch_lifts",attempted,"singular_matches",singular,
          "all_four_reconstructed",reconstructed,
          "exact_open_hits",exact_hits)
    print("Z49_STRUCTURAL_SLOPE_HENSEL_DONE")


if __name__=="__main__":
    main()
