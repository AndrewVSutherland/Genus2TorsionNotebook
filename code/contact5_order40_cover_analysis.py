#!/usr/bin/env sage --python
"""
Analysis of the remaining order-40 cover in the contact-5/order-20 family.

Usage:

    sage --python code/contact5_order40_cover_analysis.py summary
    sage --python code/contact5_order40_cover_analysis.py equation
    sage --python code/contact5_order40_cover_analysis.py search 100

The quotient equation is obtained from the condition that z*(z+1) is a square
modulo

    4*z^4 + 4*(2-s)*z^3 + 4*(2-s)*z^2 + (4-s)*z + 1.

The eliminated D,s equation is even in D; this script works with Y = D^2.
Rational order-40 specializations require s = r^2 - 1 and Y a rational square.
"""

import sys
from math import gcd as ZZ_gcd

from sage.all import GF, PolynomialRing, QQ, ZZ, factor, gcd, vector


def quotient_polynomial():
    coeff_ring = PolynomialRing(QQ, ["s", "A", "B", "C", "D"], order="lex")
    s, A, B, C, D = coeff_ring.gens()
    z_ring = PolynomialRing(coeff_ring, "z")
    z = z_ring.gen()

    quartic = (
        4 * z**4
        + 4 * (2 - s) * z**3
        + 4 * (2 - s) * z**2
        + (4 - s) * z
        + 1
    )
    target = z * (z + 1)
    square_root = A * z**3 + B * z**2 + C * z + D
    remainder = (square_root**2 - target).quo_rem(quartic)[1]
    equations = [coeff_ring(remainder[i]) for i in range(4)]

    gb_ring = PolynomialRing(QQ, ["A", "B", "C", "D", "s"], order="lex")
    A2, B2, C2, D2, s2 = gb_ring.gens()
    convert = coeff_ring.hom([s2, A2, B2, C2, D2], gb_ring)
    ideal = gb_ring.ideal([convert(e) for e in equations])
    basis = ideal.groebner_basis()
    ds_polys = [g for g in basis if set(g.variables()).issubset({D2, s2})]
    if len(ds_polys) != 1:
        raise RuntimeError("Expected a unique eliminated D,s equation.")

    ds_poly = ds_polys[0]
    quotient_ring = PolynomialRing(QQ, ["Y", "s"], order="lex")
    Y, sy = quotient_ring.gens()
    quotient = quotient_ring(0)
    for mon, coeff in ds_poly.dict().items():
        if mon[3] % 2 != 0:
            raise RuntimeError("Eliminated equation is not even in D.")
        quotient += QQ(coeff) * Y ** (mon[3] // 2) * sy ** mon[4]
    return quotient


def rational_parameters_of_height(bound):
    vals = []
    seen = set()
    for den in range(1, bound + 1):
        for num in range(-bound, bound + 1):
            if ZZ_gcd(num, den) != 1:
                continue
            r = QQ(num) / QQ(den)
            if r not in seen:
                seen.add(r)
                vals.append(r)
    return vals


def is_square_q(q):
    q = QQ(q)
    if q < 0:
        return False
    return ZZ(q.numerator()).is_square() and ZZ(q.denominator()).is_square()


def sqrt_q(q):
    q = QQ(q)
    return QQ(ZZ(q.numerator()).sqrt()) / QQ(ZZ(q.denominator()).sqrt())


def newton_summary(poly):
    pts = [tuple(mon) for mon in poly.dict()]
    vertices = [(0, 0), (2, 0), (8, 6), (8, 18), (0, 10)]
    area2 = 0
    boundary = 0
    for p1, p2 in zip(vertices, vertices[1:] + vertices[:1]):
        x1, y1 = p1
        x2, y2 = p2
        area2 += x1 * y2 - x2 * y1
        boundary += ZZ_gcd(abs(x2 - x1), abs(y2 - y1))
    interior = (abs(area2) - boundary + 2) // 2
    return vertices, boundary, interior


def print_summary(poly):
    Y, s = poly.parent().gens()
    print("Quotient equation Q(Y,s)=0 with Y=D^2")
    print("  deg_Y", poly.degree(Y), "deg_s", poly.degree(s),
          "total_degree", poly.total_degree(), "terms", len(poly.monomials()))
    print("  irreducible_factors", len(poly.factor()))
    print("")

    for sval in [0, 3, -1]:
        univar = PolynomialRing(QQ, "Y")
        Yu = univar.gen()
        spec = univar(poly.subs({s: QQ(sval)}))
        print("Special fiber s =", sval)
        print("  degree", spec.degree(), "rational_roots", spec.roots(QQ))
        print("  factor", spec.factor())
    print("")

    disc = poly.discriminant(Y)
    print("Discriminant in Y:")
    print("  degree", disc.degree(s), "terms", len(disc.monomials()))
    for f, e in disc.factor():
        label = f if f.degree(s) <= 20 else "large"
        print("  factor exp", e, "degree", f.degree(s), "terms", len(f.monomials()), label)
    print("")

    sing_ideal = poly.parent().ideal([poly, poly.derivative(Y), poly.derivative(s)])
    print("Affine singular ideal:")
    print("  dimension", sing_ideal.dimension(),
          "degree", sing_ideal.vector_space_dimension(),
          "rational_points", sing_ideal.variety(QQ))
    s_elim = [g for g in sing_ideal.groebner_basis() if set(g.variables()).issubset({s})]
    if s_elim:
        print("  s-projection:")
        for f, e in s_elim[0].factor():
            print("    exp", e, "degree", f.degree(s), "terms", len(f.monomials()))
    print("")

    vertices, boundary, interior = newton_summary(poly)
    print("Newton polygon:")
    print("  vertices", vertices)
    print("  boundary_lattice_points", boundary, "interior_lattice_points", interior)
    print("  expected toric genus before affine nodes:", interior)
    print("")

    edges = [((0, 0), (2, 0)), ((2, 0), (8, 6)), ((8, 6), (8, 18)),
             ((8, 18), (0, 10)), ((0, 10), (0, 0))]
    T_ring = PolynomialRing(QQ, "T")
    T = T_ring.gen()
    print("Newton edge checks:")
    for idx, (p1, p2) in enumerate(edges):
        x1, y1 = p1
        x2, y2 = p2
        dx = x2 - x1
        dy = y2 - y1
        length = ZZ_gcd(abs(dx), abs(dy))
        stepx = dx // length
        stepy = dy // length
        coeffs = []
        for k in range(length + 1):
            mon = (x1 + k * stepx, y1 + k * stepy)
            coeffs.append(poly.monomial_coefficient(Y ** mon[0] * s ** mon[1]))
        edge_poly = sum(QQ(c) * T**i for i, c in enumerate(coeffs))
        print("  edge", idx, p1, p2, "length", length,
              "squarefree", gcd(edge_poly, edge_poly.derivative()).degree() == 0)
        print("    ", edge_poly.factor())


def search(poly, bound):
    Y, s = poly.parent().gens()
    univar = PolynomialRing(QQ, "Y")
    params = rational_parameters_of_height(bound)
    checked = 0
    fibers_with_roots = 0
    square_hits = []
    print("Searching r of height <=", bound, "parameters", len(params))
    for r in params:
        if r in {QQ(-1), QQ(0), QQ(1)}:
            continue
        sval = r**2 - 1
        if sval in {QQ(0), QQ(-1)}:
            continue
        checked += 1
        spec = univar(poly.subs({s: sval}))
        roots = spec.roots(QQ)
        if not roots:
            continue
        fibers_with_roots += 1
        for yval, mult in roots:
            if not is_square_q(yval):
                continue
            dval = sqrt_q(yval)
            tval = (3 - r**2) / (r**2 - 1)
            square_hits.append((r, sval, yval, dval, tval, mult))
            print("HIT", "r", r, "s", sval, "Y", yval,
                  "D", dval, "t", tval, "mult", mult)
    print("DONE search height", bound)
    print("checked", checked, "fibers_with_rational_Y_roots", fibers_with_roots,
          "square_hits", len(square_hits))


def finite_sieve(poly, primes):
    Y, s = poly.parent().gens()
    print("Finite residue counts for Q(Y,s)=0 with s+1 and Y squares")
    for p in primes:
        F = GF(p)
        count = 0
        good_s = set()
        for sv in F:
            if sv in {F(0), F(-1)}:
                continue
            if not (sv + 1).is_square():
                continue
            for yv in F:
                if yv == 0 or not yv.is_square():
                    continue
                val = F(0)
                for mon, coeff in poly.dict().items():
                    val += F(coeff) * (yv ** mon[0]) * (sv ** mon[1])
                if val == 0:
                    count += 1
                    good_s.add(int(sv))
        print("p", p, "points", count, "s_residues", len(good_s))


def main():
    command = sys.argv[1] if len(sys.argv) >= 2 else "summary"
    poly = quotient_polynomial()
    if command == "summary":
        print_summary(poly)
    elif command == "equation":
        print(poly)
    elif command == "search":
        bound = int(sys.argv[2]) if len(sys.argv) >= 3 else 100
        search(poly, bound)
    elif command == "finite":
        finite_sieve(poly, [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43])
    else:
        raise SystemExit("unknown command: %s" % command)


if __name__ == "__main__":
    main()
