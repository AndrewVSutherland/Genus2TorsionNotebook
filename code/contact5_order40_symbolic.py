#!/usr/bin/env sage --python
"""
Symbolic descent setup for the contact-5/order-20 family.

Run from torsion_jac with

    sage --python code/contact5_order40_symbolic.py

The contact-5 family has

    h = 1 + t*x + ((t^2 - 1)/2)*x^2,
    f = h^2 - ((t + 1)^4/4)*x^5,
    H = [x^2 + 2*x/(t+1), (t+2)*x + 1].

To decide whether H is divisible by 2, use the usual 2-descent criterion:
the Mumford u-polynomial of H must be a square in the etale algebra cut out
by the branch polynomial.

The rational branch point x=1 first gives

    (t+3)/(t+1) = r^2,       t = (3-r^2)/(r^2-1).

Putting s = r^2 - 1 and scaling x = s*z, the quartic part becomes

    4*z^4 + 4*(2-s)*z^3 + 4*(2-s)*z^2 + (4-s)*z + 1,

and the remaining square condition is that z*(z+1) be a square modulo this
quartic.
"""

from sage.all import PolynomialRing, QQ


def scaled_quartic_data():
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
    return coeff_ring, z, quartic, target


def degree_two_component():
    coeff_ring = PolynomialRing(QQ, ["s", "B", "C", "D"], order="lex")
    s, B, C, D = coeff_ring.gens()
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
    square_root = B * z**2 + C * z + D
    remainder = (square_root**2 - target).quo_rem(quartic)[1]
    equations = [coeff_ring(remainder[i]) for i in range(4)]

    gb_ring = PolynomialRing(QQ, ["B", "C", "D", "s"], order="lex")
    B2, C2, D2, s2 = gb_ring.gens()
    convert = coeff_ring.hom([s2, B2, C2, D2], gb_ring)
    ideal = gb_ring.ideal([convert(e) for e in equations])
    return ideal.groebner_basis()


def full_cover_summary():
    coeff_ring, z, quartic, target = scaled_quartic_data()
    s, A, B, C, D = coeff_ring.gens()
    square_root = A * z**3 + B * z**2 + C * z + D
    remainder = (square_root**2 - target).quo_rem(quartic)[1]
    equations = [coeff_ring(remainder[i]) for i in range(4)]

    gb_ring = PolynomialRing(QQ, ["A", "B", "C", "D", "s"], order="lex")
    A2, B2, C2, D2, s2 = gb_ring.gens()
    convert = coeff_ring.hom([s2, A2, B2, C2, D2], gb_ring)
    ideal = gb_ring.ideal([convert(e) for e in equations])
    basis = ideal.groebner_basis()

    ds_polys = [g for g in basis if set(g.variables()).issubset({D2, s2})]
    return basis, ds_polys


print("First square condition:")
print("  (t+3)/(t+1) = r^2, so t = (3-r^2)/(r^2-1)")
print("  Put s = r^2 - 1 and x = s*z.")
print("")

print("Degree-2 square-root component, S = B*z^2 + C*z + D:")
for g in degree_two_component():
    print(" ", g.factor())
print("")
print("Consequences:")
print("  (s-3)*s*(s+1)=0")
print("  s=3 gives r^2=4 and t=-1/3.")
print("  s=0 is the t=infinity boundary; s=-1 gives t=-3, a singular member.")
print("  At s=3, one square root is S=-2*z^2 + z + 1.")
print("")

basis, ds_polys = full_cover_summary()
print("Full square-root cover:")
print("  Groebner basis length:", len(basis))
print("  polynomials involving only D and s:", len(ds_polys))
for g in ds_polys:
    print(
        "  D/s plane equation: total_degree=%s deg_D=%s deg_s=%s terms=%s factors=%s"
        % (
            g.total_degree(),
            g.degree(basis[0].parent().gen(3)),
            g.degree(basis[0].parent().gen(4)),
            len(g.monomials()),
            len(g.factor()),
        )
    )
