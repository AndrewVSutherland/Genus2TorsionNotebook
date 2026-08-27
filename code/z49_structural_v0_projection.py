#!/usr/bin/env python3
"""Exact low-memory elimination of the constant-B (v=0) order-49 slice.

Put w=u^2 (open condition w!=0), solve the second reduced equation for a,
and call the remaining equations N0,N2,N3 in Q[b,w,r].  Two w-resultants
split into a boundary, a degree-drop conic, and factors F8,G10.  Eliminating
b between F8,G10 and separately resolving the conic proves that this open
slice has no Q-point.

The process is capped below 512 MB and does not compute a Groebner basis.
"""

import argparse
import resource
import sys

import sympy as sp

sys.path.insert(0, "code")
import z49_structural_contact_iterate as contact


def derive_v0():
    _x, variables, _h, _f, _solved, equations = contact.derive_system()
    a,b,u,v,r = variables
    w = sp.symbols("w")
    ev = [poly.as_expr().subs(v,0) for poly in equations]
    asol = sp.solve(ev[1],a)[0]
    nn = []
    for index in (0,2,3):
        numerator = sp.together(ev[index].subs(a,asol)).as_numer_denom()[0]
        converted = 0
        for monomial,coefficient in sp.Poly(numerator,u,b,r).terms():
            assert monomial[0] % 2 == 0
            converted += (coefficient*w**(monomial[0]//2)*
                          b**monomial[1]*r**monomial[2])
        nn.append(sp.Poly(converted,b,w,r).primitive()[1].as_expr())
    return b,w,r,nn


def nontrivial_factor(poly, bdegree):
    for factor, exponent in sp.factor_list(poly)[1]:
        if sp.degree(factor,sp.symbols("b")) == bdegree:
            assert exponent == 1
            return factor
    raise RuntimeError("expected factor not found")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--memory-mb",type=int,default=480)
    args = parser.parse_args()
    resource.setrlimit(resource.RLIMIT_AS,
                       (args.memory_mb*1024**2,args.memory_mb*1024**2))
    b,w,r,nn = derive_v0()
    n0,n2,n3 = nn
    rw0 = sp.Poly(sp.resultant(n0,n3,w),b,r)
    rw2 = sp.Poly(sp.resultant(n2,n3,w),b,r)
    fac0 = sp.factor_list(rw0.as_expr())[1]
    fac2 = sp.factor_list(rw2.as_expr())[1]
    summary0 = [(sp.degree(f,b),sp.degree(f,r),e) for f,e in fac0]
    summary2 = [(sp.degree(f,b),sp.degree(f,r),e) for f,e in fac2]
    f8 = next(f for f,e in fac0 if sp.degree(f,b) == 8)
    g10 = next(f for f,e in fac2 if sp.degree(f,b) == 10)

    projection = sp.Poly(sp.resultant(f8,g10,b),r)
    projection_factors = sp.factor_list(projection.as_expr())[1]
    projection_degrees = sorted((sp.degree(f,r),e)
                                for f,e in projection_factors)
    assert projection.degree() == 194
    assert projection_degrees == [(1,20),(56,1),(118,1)]
    linear = next(f for f,e in projection_factors if sp.degree(f,r) == 1)
    assert sp.factor(linear) in (r-1,1-r)

    # The common conic factor in the two resultants is a leading-degree
    # artifact, so it must be substituted and checked directly.
    conic_r = (13-2*b**2)/7
    conic = [sp.cancel(poly.subs(r,conic_r)) for poly in nn]
    c03 = sp.Poly(sp.resultant(conic[0],conic[2],w),b)
    c23 = sp.Poly(sp.resultant(conic[1],conic[2],w),b)
    cf03 = sp.factor_list(c03.as_expr())[1]
    cf23 = sp.factor_list(c23.as_expr())[1]
    deg03 = sorted((sp.degree(f,b),e) for f,e in cf03)
    deg23 = sorted((sp.degree(f,b),e) for f,e in cf23)
    assert deg03 == [(2,5),(15,1)]
    assert deg23 == [(2,10),(25,1)]
    p15 = next(f for f,e in cf03 if sp.degree(f,b) == 15)
    p25 = next(f for f,e in cf23 if sp.degree(f,b) == 25)
    assert sp.gcd(sp.Poly(p15,b),sp.Poly(p25,b)).degree() == 0
    quadratic = next(f for f,e in cf03 if sp.degree(f,b) == 2)
    assert sp.factor(quadratic) in (b**2-3,3-b**2)
    failed_denominator = sp.Poly(conic[2],w).LC()
    denominator_factors = sp.factor_list(failed_denominator)[1]
    assert len(denominator_factors) == 1
    assert sp.degree(denominator_factors[0][0],b) == 9

    print("Z49_STRUCTURAL_V0_PROJECTION")
    print("w_resultant_factors_0",summary0)
    print("w_resultant_factors_2",summary2)
    print("main_projection_degree",projection.degree(),
          "factor_degrees",projection_degrees)
    print("only_rational_linear_projection_root r=1 excluded_open true")
    print("conic_backsub_factor_degrees",deg03,deg23,
          "P15_P25_gcd_degree",0)
    print("conic_rational_quadratic b^2-3 irreducible true")
    print("failed_denominator_degree 9 irreducible true")
    print("CONCLUSION open_v0_Q_points 0")
    print("Z49_STRUCTURAL_V0_PROJECTION_DONE")


if __name__ == "__main__":
    main()
