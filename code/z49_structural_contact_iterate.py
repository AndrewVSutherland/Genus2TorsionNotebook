#!/usr/bin/env python3
"""Derive and scout the iterated contact-7 equations for order 49.

Start with

    h = 1-7*x/2+a*x^2+b*x^3,
    f = (h^2+(x-1)^7)/x^2,

so P=(1,h(1)) gives D=P-infinity of order 7.  If

    psi = (A+B*y)/(x-1),  deg(A)=4, deg(B)<=1,

has divisor 7*R-P-6*infinity, then 7*(R-infinity)=D.  Its norm identity is

    A^2-B^2*f = (x-1)*(x-r)^7.

After making A monic, coefficients x^7,x^6,x^5,x^4 solve successively for
the four remaining coefficients of A.  This leaves four explicit equations
in (a,b,B0,B1,r).  The script derives those equations exactly, reports their
sizes, performs small finite-field open/smooth diagnostics, and checks a
small integral box.  It uses no Groebner basis.
"""

from itertools import product
import argparse

import sympy as sp


def primitive_integer_poly(expr, variables):
    poly = sp.Poly(sp.cancel(expr), *variables, domain=sp.QQ)
    den = sp.ilcm(*[term.q for term in poly.coeffs()]) if poly.coeffs() else 1
    out = sp.Poly(poly.as_expr() * den, *variables, domain=sp.ZZ)
    content, primitive = sp.polys.polytools.primitive(out)
    if primitive.LC() < 0:
        primitive = -primitive
    return primitive


def derive_system():
    x = sp.symbols("x")
    a, b, u, v, r = sp.symbols("a b u v r")
    z0, z1, z2, z3 = sp.symbols("z0 z1 z2 z3")
    h = 1 - sp.Rational(7, 2)*x + a*x**2 + b*x**3
    numerator = sp.expand(h**2 + (x-1)**7)
    assert sp.rem(numerator, x**2, domain=sp.QQ[a,b]) == 0
    f = sp.cancel(numerator/x**2)
    A = x**4 + z3*x**3 + z2*x**2 + z1*x + z0
    B = u + v*x
    error = sp.Poly(sp.expand(A**2-B**2*f-(x-1)*(x-r)**7), x)

    solved = {}
    current = error.as_expr()
    for degree, variable in zip((7, 6, 5, 4), (z3, z2, z1, z0)):
        coefficient = sp.Poly(current, x).coeff_monomial(x**degree)
        coefficient = sp.cancel(coefficient.subs(solved))
        answer = sp.solve(coefficient, variable, dict=False)
        assert len(answer) == 1
        solved[variable] = sp.factor(answer[0])
        current = sp.cancel(current.subs(variable, solved[variable]))

    variables = (a, b, u, v, r)
    residuals = []
    pcurrent = sp.Poly(sp.together(current), x)
    for degree in range(4):
        residuals.append(primitive_integer_poly(
            pcurrent.coeff_monomial(x**degree), variables))
    assert all(sp.cancel(sp.Poly(current, x).coeff_monomial(x**i)) == 0
               for i in range(4, 9))
    return x, variables, h, f, solved, residuals


def term_data(poly):
    return [(mon, int(coeff)) for mon, coeff in poly.terms()]


def eval_terms(terms, values, p=None):
    out = 0
    for monomial, coeff in terms:
        z = coeff
        for value, exponent in zip(values, monomial):
            z *= value**exponent
        out += z
    return out if p is None else out % p


def trim_mod(poly, p):
    poly = [c % p for c in poly]
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def gcd_mod(a, b, p):
    a, b = trim_mod(a, p), trim_mod(b, p)
    while b != [0]:
        rr = a[:]
        while len(rr) >= len(b) and rr != [0]:
            shift = len(rr)-len(b)
            scale = rr[-1]*pow(b[-1], -1, p) % p
            for i in range(len(b)):
                rr[i+shift] = (rr[i+shift]-scale*b[i]) % p
            rr = trim_mod(rr, p)
        a, b = b, rr
    if a == [0]:
        return a
    inv = pow(a[-1], -1, p)
    return [(c*inv) % p for c in a]


def f_coefficients_mod(a, b, p):
    # Derived directly from (h^2+(x-1)^7)/x^2, low to high.
    inv2 = pow(2, -1, p)
    h = [1, (-7*inv2) % p, a % p, b % p]
    hs = [0]*7
    for i, hi in enumerate(h):
        for j, hj in enumerate(h):
            hs[i+j] = (hs[i+j]+hi*hj) % p
    seventh = [sp.binomial(7, i)*((-1)**(7-i)) % p for i in range(8)]
    num = [0]*8
    for i in range(8):
        num[i] = ((hs[i] if i < 7 else 0)+int(seventh[i])) % p
    assert num[0] == num[1] == 0
    return trim_mod(num[2:], p)


def smooth_mod(a, b, p):
    f = f_coefficients_mod(a, b, p)
    derivative = [(i*f[i]) % p for i in range(1, len(f))]
    return len(f) == 6 and len(gcd_mod(f, derivative, p)) == 1


def finite_scout(polys, primes):
    terms = [term_data(poly) for poly in polys]
    for p in primes:
        if p == 2:
            continue
        raw = opened = smooth = 0
        samples = []
        inv2 = pow(2, -1, p)
        for values in product(range(p), repeat=5):
            if any(eval_terms(eq, values, p) for eq in terms):
                continue
            raw += 1
            a, b, u, v, r = values
            hp = (a+b-5*inv2) % p
            if r == 1 or hp == 0 or (u+v) % p == 0 or (u+v*r) % p == 0:
                continue
            opened += 1
            if not smooth_mod(a, b, p):
                continue
            smooth += 1
            if len(samples) < 4:
                samples.append(values)
        print("finite", p, "raw", raw, "open", opened,
              "open_smooth", smooth, "samples", samples)


def integral_scout(polys, height):
    terms = [term_data(poly) for poly in polys]
    checked = raw = opened = 0
    samples = []
    for values in product(range(-height, height+1), repeat=5):
        checked += 1
        if any(eval_terms(eq, values) for eq in terms):
            continue
        raw += 1
        a, b, u, v, r = values
        if r == 1 or 2*(a+b)-5 == 0 or u+v == 0 or u+v*r == 0:
            continue
        opened += 1
        if len(samples) < 10:
            samples.append(values)
    print("integral_box", height, "checked", checked, "raw", raw,
          "open", opened, "samples", samples)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=3)
    parser.add_argument("--primes", default="3,5,11")
    parser.add_argument("--print-equations", action="store_true")
    args = parser.parse_args()
    _x, variables, _h, _f, solved, residuals = derive_system()
    print("Z49_STRUCTURAL_CONTACT_ITERATE")
    print("variables", variables)
    for variable in (sp.symbols("z3"), sp.symbols("z2"),
                     sp.symbols("z1"), sp.symbols("z0")):
        formula = solved[variable]
        print("triangular", variable, "total_degree",
              sp.Poly(formula.as_numer_denom()[0], *variables).total_degree(),
              "terms", len(sp.Poly(formula.as_numer_denom()[0], *variables).terms()),
              "denominator", formula.as_numer_denom()[1])
        if args.print_equations:
            print(variable, "=", formula)
    for i, poly in enumerate(residuals):
        print("equation", i, "total_degree", poly.total_degree(),
              "terms", len(poly.terms()))
        if args.print_equations:
            print("E%d =" % i, sp.factor(poly.as_expr()))
    finite_scout(residuals,
                 [int(q) for q in args.primes.split(",") if q.strip()])
    integral_scout(residuals, args.height)
    print("Z49_STRUCTURAL_CONTACT_ITERATE_DONE")


if __name__ == "__main__":
    main()
