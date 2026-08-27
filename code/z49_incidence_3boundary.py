#!/usr/bin/env python3
"""Resolve and search the first 3-adic boundary charts of the [49] incidence.

The four incidence equations are derived by
``z49_structural_contact_iterate.py``.  There are two open points modulo 3,

    (a,b,u,v,r) = (1,1,eps,0,0),  eps = +/-1.

This script does four bounded-memory tasks.

* It derives the tangent space and the first weighted exceptional divisor at
  both points.
* It follows the pullback of the Q5 discriminant component through several
  exact modular Hensel layers.
* It searches the resulting one-dimensional branches using the local
  parameter T=(b-1)/3, rather than a five-variable box.
* It derives three initial-form stages on the balanced first-pole chart
  v3(b)=v3(v)=-1 with a,u,r integral.  That chart is still nonreduced after
  the displayed stages, so the script deliberately does not pretend to
  search it exhaustively.

Every reported rational hit is substituted into all four original equations
and all open conditions.  The default run stays well below 200 MB RSS.
"""

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import product
from math import gcd
import argparse
import sys

import sympy as sp

sys.path.insert(0, "code")
import z49_structural_contact_iterate as contact
import z49_structural_3adic as hensel


def q5(a, b):
    """The non-contact component of 256*Disc(f)=H^7*Q5."""
    return (
        432*a**4 - 64*a**3*b**3 + 1008*a**3*b + 3024*a**3
        - 448*a**2*b**3 + 224*a**2*b**2 + 21168*a**2*b
        - 32536*a**2 - 2016*a*b**4 + 4480*a*b**3
        + 38416*a*b**2 - 109760*a*b + 78890*a - 864*b**5
        + 5936*b**4 + 7056*b**3 - 96040*b**2 + 120050*b - 60025
    )


def v3_integer(value):
    if value == 0:
        return 10**9
    value = abs(int(value))
    answer = 0
    while value % 3 == 0:
        value //= 3
        answer += 1
    return answer


def v3_fraction(value):
    value = Fraction(value)
    if value == 0:
        return 10**9
    return v3_integer(value.numerator)-v3_integer(value.denominator)


def residue_fraction(value, modulus):
    value = Fraction(value)
    return value.numerator*pow(value.denominator, -1, modulus) % modulus


def signed_mod3(value):
    value %= 3
    return -1 if value == 2 else value


def exceptional_equations(equations, variables, eps):
    """Return E_i/9 mod 3 in the weighted chart, without large expansion."""
    a, b, u, v, r = variables
    A, T, U, V, R = sp.symbols("A T U V R")
    base = {a: 1, b: 1, u: eps, v: 0, r: 0}
    answer = []
    for equation in equations:
        constant = int(equation.eval(base)) // 9
        linear_b = int(equation.diff(b).eval(base)) // 3
        # Coefficient of (b-1)^2, avoiding division by a symbolic factorial.
        z = sp.symbols("z")
        in_b = sp.Poly(equation.as_expr().subs({a: 1, b: 1+z,
                                                u: eps, v: 0, r: 0}), z)
        quadratic_b = int(in_b.coeff_monomial(z**2))
        expr = constant + linear_b*T + quadratic_b*T**2
        for variable, new_variable in ((a,A),(u,U),(v,V),(r,R)):
            expr += int(equation.diff(variable).eval(base))*new_variable
        expr = sp.Poly(expr, A,T,U,V,R, modulus=3).as_expr()
        answer.append(expr)
    return (A,T,U,V,R), answer


def derive_integral_chart(equations, variables):
    a, b, u, v, r = variables
    unknown_indices = [0,2,3,4]
    branches = []
    print("INTEGRAL_OPEN_CHARTS")
    for eps in (1,-1):
        base = {a:1, b:1, u:eps % 3, v:0, r:0}
        jacobian = [[int(eq.diff(variables[j]).eval(base)) % 3
                     for j in range(5)] for eq in equations]
        transverse = [[row[j] for j in unknown_indices] for row in jacobian]
        inverse = hensel.matrix_inverse_mod(transverse, 3)
        assert inverse is not None
        # db=1 and all other tangent coordinates zero.
        tangent_rhs = [-row[1] % 3 for row in jacobian]
        assert hensel.mat_vec(inverse, tangent_rhs, 3) == [0,0,0,0]
        print("eps",eps,"jacobian_mod3",jacobian)
        print("eps",eps,"tangent_db_1",[0,1,0,0,0])

        new_variables, exceptional = exceptional_equations(
            equations, variables, eps)
        A,T,U,V,R = new_variables
        expected = [
            A+T**2+T+1,
            A+T**2+eps*U+eps*V-1,
            -A+R-T**2-eps*U+eps*V+1,
            A+R+T**2-T-eps*U-eps*V+1,
        ]
        assert all(sp.Poly(g-h,A,T,U,V,R,modulus=3).is_zero
                   for g,h in zip(exceptional,expected))
        line = {A:-(T**2+T+1), U:eps*T, V:-eps, R:-1}
        assert all(sp.Poly(g.subs(line),T,modulus=3).is_zero
                   for g in exceptional)
        # There are 3^5 points and only the three points on this line.
        finite_solutions = []
        for values in product(range(3), repeat=5):
            sub = dict(zip(new_variables, values))
            if all(int(g.subs(sub)) % 3 == 0 for g in exceptional):
                finite_solutions.append(values)
        assert finite_solutions == sorted([
            ((-(t*t+t+1)) % 3,t,(eps*t) % 3,(-eps) % 3,2)
            for t in range(3)])
        print("eps",eps,"exceptional_equations",[str(g) for g in expected])
        print("eps",eps,
              "exceptional_line A=-(T^2+T+1) U=eps*T V=-eps R=-1")
        branches.append((eps,inverse))
    return branches


def modular_q5_layers(polys, maxima, branches, layers):
    """Certify successive residue equations for Q5 on both Hensel branches."""
    print("Q5_EXCEPTIONAL_LAYERS")
    center = 0
    scale = 1
    for level in range(layers):
        depth = 7*(level+1)
        branch_rows = []
        for eps, inverse in branches:
            quotients = []
            for digit in range(3):
                T = center+scale*digit
                bvalue = 1+3*T
                residues, modulus = hensel.lift(
                    polys,maxima,inverse,eps,(bvalue,1),depth+1)
                qvalue = q5(residues[0],bvalue) % modulus
                assert qvalue % (3**depth) == 0
                quotients.append((qvalue//(3**depth)) % 3)
            branch_rows.append(quotients)
        assert branch_rows[0] == branch_rows[1]
        quotients = branch_rows[0]
        roots = [digit for digit,value in enumerate(quotients) if value == 0]
        assert len(roots) == 1
        # Fit q(W)=constant+slope*W over F3 for a compact equation.
        constant = quotients[0]
        slope = (quotients[1]-constant) % 3
        assert all((constant+slope*w) % 3 == quotients[w]
                   for w in range(3))
        print("level",level+1,"T=center+scale*W",
              "center",center,"scale",scale,"depth",depth,
              "quotients",quotients,
              "equation",f"{signed_mod3(constant)}"
              f"{signed_mod3(slope):+d}*W",
              "deeper_digit",roots[0])
        center += scale*roots[0]
        scale *= 3
    print("common_deeper_disc T=",center,"mod",scale)


def rational_parameters(height):
    """Reduced 3-integral rationals of local-coordinate height <= height."""
    for denominator in range(1,height+1):
        if denominator % 3 == 0:
            continue
        for numerator in range(-height,height+1):
            if gcd(abs(numerator),denominator) == 1:
                yield numerator,denominator


def exact_open(values, poly_terms):
    a,b,u,v,r = values
    if not all(hensel.evaluate_exact(poly,values) == 0 for poly in poly_terms):
        return False
    return (r != 1 and a+b != Fraction(5,2) and u+v != 0
            and u+v*r != 0 and q5(a,b) != 0)


def local_parameter_search(polys, maxima, branches, height, precision,
                           coordinate_bound):
    assert 3**precision > 2*coordinate_bound**2
    tvalues = list(rational_parameters(height))
    reconstructed = 0
    exact_hits = []
    thickness = Counter()
    for tfrac in tvalues:
        tq = Fraction(*tfrac)
        bq = 1+3*tq
        bfrac = (bq.numerator,bq.denominator)
        for eps,inverse in branches:
            residues,modulus = hensel.lift(
                polys,maxima,inverse,eps,bfrac,precision)
            bresidue = residue_fraction(bq,modulus)
            qvalue = q5(residues[0],bresidue) % modulus
            valuation = v3_integer(qvalue)
            thickness[(f">={precision}" if valuation >= precision else valuation)] += 1
            rec = [hensel.rational_reconstruction_bounded(
                       value,modulus,coordinate_bound) for value in residues]
            if any(value is None for value in rec):
                continue
            reconstructed += 1
            aq,uq,vq,rq = [Fraction(*value) for value in rec]
            values = [aq,bq,uq,vq,rq]
            if exact_open(values,polys):
                exact_hits.append(values)
                print("EXACT_OPEN_HIT",values)
    print("LOCAL_T_SEARCH height",height,"precision",precision,
          "coordinate_bound",coordinate_bound)
    print("rational_T",len(tvalues),"branches",2,
          "lifts",2*len(tvalues),"thickness",dict(thickness))
    print("all_four_reconstructed",reconstructed,
          "exact_open_hits",len(exact_hits))
    print("coverage T=(b-1)/3 3-integral local_height<=",height,
          "other_coordinate_bound<=",coordinate_bound)
    return exact_hits


def sparse_multiply(left,right):
    answer = defaultdict(Fraction)
    for lm,lc in left.items():
        for rm,rc in right.items():
            answer[tuple(x+y for x,y in zip(lm,rm))] += lc*rc
    return {m:c for m,c in answer.items() if c}


def sparse_power(poly, exponent, nvariables):
    answer = {(0,)*nvariables:Fraction(1)}
    base = poly
    while exponent:
        if exponent & 1:
            answer = sparse_multiply(answer,base)
        exponent //= 2
        if exponent:
            base = sparse_multiply(base,base)
    return answer


def initial_form(poly, substitutions, new_variables):
    """Exact sparse substitution followed by the lowest 3-adic initial form."""
    n = len(new_variables)
    maxima = [max(m[j] for m,_ in poly.terms()) for j in range(5)]
    powers = []
    for j,substitution in enumerate(substitutions):
        powers.append([sparse_power(substitution,e,n)
                       for e in range(maxima[j]+1)])
    expanded = defaultdict(Fraction)
    for monomial,coefficient in poly.terms():
        row = {(0,)*n:Fraction(int(coefficient))}
        for j,exponent in enumerate(monomial):
            row = sparse_multiply(row,powers[j][exponent])
        for new_monomial,new_coefficient in row.items():
            expanded[new_monomial] += new_coefficient
    expanded = {m:c for m,c in expanded.items() if c}
    minimum = min(v3_fraction(c) for c in expanded.values())
    expr = 0
    for monomial,coefficient in expanded.items():
        if v3_fraction(coefficient) != minimum:
            continue
        if minimum >= 0:
            unit = coefficient/Fraction(3**minimum)
        else:
            unit = coefficient*Fraction(3**(-minimum))
        residue = unit.numerator*pow(unit.denominator,-1,3) % 3
        expr += signed_mod3(residue)*sp.prod(
            variable**exponent
            for variable,exponent in zip(new_variables,monomial))
    return minimum,sp.Poly(expr,*new_variables,modulus=3).as_expr()


def pole_chart(equations):
    """Three controlled initial-form stages on one balanced first-pole chart."""
    print("BALANCED_FIRST_POLE_CHART")
    A,B,U,V,R = sp.symbols("A B U V R")
    new = (A,B,U,V,R)
    zero = (0,0,0,0,0)
    unit = lambda index,coefficient=1: {
        tuple(1 if i == index else 0 for i in range(5)):Fraction(coefficient)}
    chart1 = [unit(0),unit(1,Fraction(1,3)),unit(2),
              unit(3,Fraction(1,3)),unit(4)]
    first = [initial_form(eq,chart1,new) for eq in equations]
    print("stage1 substitutions a=A b=B/3 u=U v=V/3 r=R; B,V units")
    print("stage1",[(valuation,str(expr)) for valuation,expr in first])
    # E3 is V^8*(V^2-B^2), so V=delta*B over F3.
    for delta in (1,-1):
        A,B,U,W,R = sp.symbols("A B U W R")
        new2 = (A,B,U,W,R)
        chart2 = [unit(0),unit(1,Fraction(1,3)),unit(2),
                  { (0,1,0,0,0):Fraction(delta,3),
                    (0,0,0,1,0):Fraction(1)},unit(4)]
        second = [initial_form(eq,chart2,new2) for eq in equations]
        print("stage2 delta",delta,"v=delta*B/3+W")
        print("stage2",[(valuation,str(expr)) for valuation,expr in second])
        # The final equation forces W=-delta*B.  Blow it up once more.
        A,B,U,X,R = sp.symbols("A B U X R")
        new3 = (A,B,U,X,R)
        chart3 = [unit(0),unit(1,Fraction(1,3)),unit(2),
                  { (0,1,0,0,0):Fraction(-2*delta,3),
                    (0,0,0,1,0):Fraction(3)},unit(4)]
        third = [initial_form(eq,chart3,new3) for eq in equations]
        print("stage3 delta",delta,
              "v=delta*B/3-delta*B+3*X")
        print("stage3",[(valuation,str(expr)) for valuation,expr in third])
        divisor = sp.Poly(delta*B*X-R+1,*new3,modulus=3)
        assert all(sp.div(sp.Poly(expr,*new3,modulus=3),divisor)[1].is_zero
                   for _,expr in third)
        print("stage3_reduced_relation",f"{delta}*B*X-R+1=0",
              "A,U,R remain free on this initial divisor")
    print("pole_scope v3(b)=v3(v)=-1; v3(a),v3(u),v3(r)>=0 only")
    print("pole_limitation chart remains nonreduced/not low-dimensional; no box search")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height",type=int,default=100)
    parser.add_argument("--precision",type=int,default=18)
    parser.add_argument("--coordinate-bound",type=int,default=5000)
    parser.add_argument("--layers",type=int,default=6)
    parser.add_argument("--skip-search",action="store_true")
    args = parser.parse_args()

    _x,variables,_h,_f,_solved,equations = contact.derive_system()
    polys = [hensel.term_data(poly) for poly in equations]
    maxima = [max(monomial[j] for terms in polys for monomial,_ in terms)
              for j in range(5)]
    print("Z49_INCIDENCE_3BOUNDARY")
    branches = derive_integral_chart(equations,variables)
    modular_q5_layers(polys,maxima,branches,args.layers)
    if not args.skip_search:
        local_parameter_search(polys,maxima,branches,args.height,
                               args.precision,args.coordinate_bound)
    pole_chart(equations)
    print("Z49_INCIDENCE_3BOUNDARY_DONE")


if __name__ == "__main__":
    main()
