#!/usr/bin/env python3
"""Root-centred normalization of the nonconstant-B order-49 incidence.

Write B=v*(x-q).  Evaluating the norm identity at x=q shows that
(q-r)/(q-1)=t^2 on the open chart, so r=q-(q-1)*t^2.  In X=x-q,
the constant and linear coefficients determine A(q),A'(q); the next two
coefficients then determine the remaining lower Taylor coefficients of A.
Only the X^4,...,X^7 equations remain.  This avoids a Groebner basis.
"""

import argparse
from fractions import Fraction
import resource
import sympy as sp


def primitive(expr,variables):
    poly=sp.Poly(sp.cancel(expr),*variables,domain=sp.QQ)
    denominator=sp.ilcm(*[coefficient.q for coefficient in poly.coeffs()])
    poly=sp.Poly(poly.as_expr()*denominator,*variables,domain=sp.ZZ)
    _content,poly=sp.polys.polytools.primitive(poly)
    return -poly if poly.LC()<0 else poly


def derive(qvalue,sign=1):
    X=sp.symbols("X");a,b,w,t=sp.symbols("a b w t")
    q=sp.Rational(qvalue.numerator,qvalue.denominator);x=X+q
    h=1-sp.Rational(7,2)*x+a*x**2+b*x**3
    f=sp.Poly(sp.cancel((h**2+(x-1)**7)/x**2),X)
    r=q-(q-1)*t**2
    rhs=sp.Poly(sp.expand((x-1)*(x-r)**7),X)
    fc=[f.coeff_monomial(X**i) for i in range(6)]
    rc=[rhs.coeff_monomial(X**i) for i in range(9)]
    c0=sign*(q-1)**4*t**7
    c1=sign*(q-1)**3*t**5*(t**2+7)/2
    c2=sp.cancel((w*fc[0]+rc[2]-c1**2)/(2*c0))
    c3=sp.cancel((w*fc[1]+rc[3]-2*c1*c2)/(2*c0))
    residuals=[
        c2**2+2*c1*c3+2*c0-w*fc[2]-rc[4],
        2*c2*c3+2*c1-w*fc[3]-rc[5],
        c3**2+2*c2-w*fc[4]-rc[6],
        2*c3-w*fc[5]-rc[7],
    ]
    variables=(a,b,w,t)
    equations=[primitive(sp.together(expr).as_numer_denom()[0],variables)
               for expr in residuals]
    A=sp.expand(X**4+c3*X**3+c2*X**2+c1*X+c0)
    error=sp.Poly(sp.cancel(A**2-w*X**2*f.as_expr()-
                            (X+q-1)*(X+q-r)**7),X)
    assert all(sp.cancel(error.coeff_monomial(X**i))==0
               for i in (0,1,2,3,8))
    return X,variables,q,r,h,f.as_expr(),A,equations


def main():
    parser=argparse.ArgumentParser();parser.add_argument("--q",default="-1")
    parser.add_argument("--sign",type=int,choices=(-1,1),default=1)
    parser.add_argument("--memory-mb",type=int,default=300)
    parser.add_argument("--print-equations",action="store_true")
    args=parser.parse_args();resource.setrlimit(resource.RLIMIT_AS,
        (args.memory_mb*1024**2,args.memory_mb*1024**2))
    _X,variables,q,r,_h,_f,_A,equations=derive(Fraction(args.q),args.sign)
    print("Z49_NONCONSTANT_ROOT_NORMALIZATION","q",q,
          "T",(-1/q if q else "infinity"),"sign",args.sign,"r",r)
    for degree,equation in enumerate(equations,4):
        print("equation_X%d"%degree,"total_degree",equation.total_degree(),
              "degrees",[equation.degree(z) for z in variables],
              "terms",len(equation.terms()))
        if args.print_equations:print("E%d ="%degree,equation.as_expr())
    print("Z49_NONCONSTANT_ROOT_NORMALIZATION_DONE")


if __name__=="__main__":main()
