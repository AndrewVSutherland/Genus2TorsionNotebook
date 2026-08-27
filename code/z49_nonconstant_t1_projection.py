#!/usr/bin/env python3
"""Staged exact projection of the genuinely nonconstant T=v/u=1 chart.

The root normalization has q=-1.  Its X^7 equation is linear in w=v^2.
On the regular chart where that coefficient is nonzero, eliminate w, strip
the open-boundary factors t=0 and t=1, and take two resultants in a.  The
small b-degree component can optionally be projected exactly to t.  No
Groebner basis is computed.
"""

import argparse
from fractions import Fraction
import resource
import sys
import sympy as sp

sys.path.insert(0,"code")
from z49_nonconstant_root_normalization import derive


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--sign",type=int,choices=(-1,1),default=1)
    parser.add_argument("--small-projection",action="store_true")
    parser.add_argument("--memory-mb",type=int,default=300)
    args=parser.parse_args();resource.setrlimit(resource.RLIMIT_AS,
        (args.memory_mb*1024**2,args.memory_mb*1024**2))
    _X,variables,_q,_r,_h,_f,_A,equations=derive(Fraction(-1),args.sign)
    a,b,w,t=variables
    e=[poly.as_expr() for poly in equations]
    assert sp.degree(e[3],w)==1
    wsolution=sp.solve(e[3],w)[0]
    denominator=sp.Poly(sp.denom(wsolution),a,b,t)
    reduced=[]
    for equation in e[:3]:
        numerator=sp.cancel(equation.subs(w,wsolution)).as_numer_denom()[0]
        factors=sp.factor_list(numerator)[1]
        core=[factor for factor,exponent in factors
              if sp.degree(factor,a)>0 or sp.degree(factor,b)>0]
        assert len(core)==1
        reduced.append(sp.Poly(core[0],a,b,t).primitive()[1].as_expr())
    print("Z49_NONCONSTANT_T1_PROJECTION","sign",args.sign)
    print("w_denominator_degrees",denominator.degree(a),
          denominator.degree(b),denominator.degree(t))
    print("reduced_profiles",[(sp.degree(p,a),sp.degree(p,b),sp.degree(p,t),
                                len(sp.Poly(p,a,b,t).terms())) for p in reduced])
    r01=sp.factor_list(sp.resultant(reduced[0],reduced[1],a))[1]
    r02=sp.factor_list(sp.resultant(reduced[0],reduced[2],a))[1]
    summary01=[(sp.degree(f,b),sp.degree(f,t),exponent) for f,exponent in r01]
    summary02=[(sp.degree(f,b),sp.degree(f,t),exponent) for f,exponent in r02]
    print("a_resultant_01_factors",summary01)
    print("a_resultant_02_factors",summary02)
    if args.small_projection:
        f2=next(f for f,e in r01 if sp.degree(f,b)==2)
        h12=next(f for f,e in r02 if sp.degree(f,b)==12)
        projection=sp.Poly(sp.resultant(f2,h12,b),t)
        factors=sp.factor_list(projection.as_expr())[1]
        print("small_projection_degree",projection.degree(),
              "terms",len(projection.terms()),
              "factor_degrees",[(sp.degree(f,t),e) for f,e in factors])
        assert [(sp.degree(f,t),e) for f,e in factors]==[(1,40),(248,1)]
        linear=next(f for f,e in factors if sp.degree(f,t)==1)
        assert sp.factor(linear) in (t,-t)
        print("small_regular_component_open_Q_points 0")
    print("scope regular_w_elimination_chart; large_b10_b12_component unresolved")
    print("Z49_NONCONSTANT_T1_PROJECTION_DONE")


if __name__=="__main__":main()
