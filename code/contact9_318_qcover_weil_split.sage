#!/usr/bin/env sage
"""Split the contact-9 3-torsion cover using the Weil pairing with 3*D9.

For a candidate line represented by

    q=x^2+u*x+v,  G^2-z*f=q^3,

write z=s^2.  The normalized function G-s*y has divisor three times the
candidate class.  Comparing it with the corresponding marked function gives
an explicit resultant condition for pairing one.  Its norm under s -> -s is
a polynomial W(u,v,z).  On the degree-40 unsigned cover, W should select the
marked line and its 12 orthogonal companions; the complementary orbit has
degree 27.
"""

import argparse
import sys

from sage.all import PolynomialRing, QQ


def family(global_base, parameter):
    if global_base:
        parameter_ring = PolynomialRing(QQ, "t")
        base = parameter_ring.fraction_field()
        t = base(parameter_ring.gen())
    else:
        parameter_ring = None
        base = QQ
        t = QQ(parameter)

    x_ring = PolynomialRing(base, "x")
    x = x_ring.gen()
    r = 1 - t**2
    if not global_base and r == 0:
        raise ValueError("t=+-1 is excluded from this presentation")
    h_at_r = 1 - base(9)/2*r + base(63)/8*r**2 - base(105)/16*r**3
    contact_a = (t**9 - h_at_r)/r**4
    marked_h = (
        1 - base(9)/2*x + base(63)/8*x**2 - base(105)/16*x**3
        + contact_a*x**4
    )
    f, remainder = (marked_h**2 + (x - 1)**9).quo_rem(x**4)
    if remainder != 0 or f.degree() != 5 or f[5] != 1:
        raise ArithmeticError("bad contact-9 quintic")
    return parameter_ring, base, t, x_ring, f


def qcover(base, f):
    ring = PolynomialRing(base, ("u", "v", "z"), order="degrevlex")
    u, v, z = ring.gens()
    coefficients = [ring(f[index]) for index in range(6)]
    a = (3*u + z)/2
    b = (z*coefficients[4] + 3*u**2 + 3*v - a**2)/2
    c = (z*coefficients[3] + u**3 + 6*u*v - 2*a*b)/2
    equations = (
        b**2 + 2*a*c - z*coefficients[2] - 3*(u**2*v + v**2),
        2*b*c - z*coefficients[1] - 3*u*v**2,
        c**2 - z*coefficients[0] - v**3,
    )
    return ring, coefficients, (a, b, c), equations


def marked_data(base, t, coefficients):
    u0 = (
        3*t**5 + base(15)/4*t**4 - 3*t**3 - base(147)/16*t**2
        - base(27)/4*t - base(27)/16
    )/(t + 1)**4
    v0 = (
        -2*t**5 - base(47)/16*t**4 + base(1)/4*t**3
        + base(7)/2*t**2 + base(11)/4*t + base(11)/16
    )/(t + 1)**4
    s0 = (
        2*(t + base(1)/2)**2
        * (t**3 + base(3)/4*t**2 - base(1)/4)
        / (t + 1)**4
    )
    z0 = s0**2
    a0 = (3*u0 + z0)/2
    b0 = (z0*coefficients[4] + 3*u0**2 + 3*v0 - a0**2)/2
    c0 = (
        z0*coefficients[3] + u0**3 + 6*u0*v0 - 2*a0*b0
    )/2
    return u0, v0, z0, s0, (a0, b0, c0)


def pairing_polynomial(ring, triangular, marked):
    u, v, z = ring.gens()
    a, b, c = triangular
    u0, v0, z0, _s0, marked_triangular = marked
    a0, b0, c0 = marked_triangular

    lambda_ring = PolynomialRing(ring, "lam")
    lam = lambda_ring.gen()
    x_ring = PolynomialRing(lambda_ring, "X")
    X = x_ring.gen()
    q = X**2 + u*X + v
    q0 = X**2 + ring(u0)*X + ring(v0)
    G = X**3 + a*X**2 + b*X + c
    G0 = X**3 + ring(a0)*X**2 + ring(b0)*X + ring(c0)
    difference = G - lam*G0

    value_on_marked = q0.resultant(difference)
    value_on_candidate_scaled = q.resultant(difference)
    pairing_equation = lam**2*value_on_marked - value_on_candidate_scaled
    relation = lam**2 - ring(z/z0)
    _quotient, remainder = pairing_equation.quo_rem(relation)
    if remainder.degree() > 1:
        raise ArithmeticError("failed to reduce the pairing equation")
    constant = ring(remainder[0])
    linear = ring(remainder[1])
    norm = constant**2 - ring(z/z0)*linear**2
    return norm, constant, linear


def analyze(global_base, parameter, check_intersections):
    parameter_ring, base, t, _x_ring, f = family(global_base, parameter)
    ring, coefficients, triangular, equations = qcover(base, f)
    u, v, z = ring.gens()
    marked = marked_data(base, t, coefficients)
    u0, v0, z0, _s0, _marked_triangular = marked

    label = "Q(t)" if global_base else f"Q at t={t}"
    print(f"BEGIN base={label}", flush=True)
    raw = ring.ideal(equations)
    boundary = ring.ideal([z, u**2 - 4*v])
    print(f"raw_dimension={raw.dimension()}", flush=True)
    residual = raw
    for iteration in (1, 2):
        residual = residual.quotient(boundary)
        print(
            f"boundary_quotient_iteration={iteration} "
            f"dimension={residual.dimension()}",
            flush=True,
        )
    if residual.dimension() != 0:
        raise ArithmeticError("double boundary quotient is not zero-dimensional")
    residual = ring.ideal(residual.groebner_basis())
    total_degree = residual.vector_space_dimension()
    print(f"total_degree={total_degree}", flush=True)
    if total_degree != 40:
        raise ArithmeticError(f"expected degree 40, got {total_degree}")

    pairing_norm, remainder_constant, remainder_linear = pairing_polynomial(
        ring, triangular, marked
    )
    pairing_norm = residual.reduce(pairing_norm)
    print(
        f"pairing_remainder_degrees="
        f"({remainder_constant.degree()},{remainder_linear.degree()})",
        flush=True,
    )
    print(
        f"pairing_norm_degree={pairing_norm.degree()} "
        f"terms={len(pairing_norm.monomials())}",
        flush=True,
    )
    print(
        f"pairing_norm_at_marked={pairing_norm(u0, v0, z0)}",
        flush=True,
    )
    if pairing_norm(u0, v0, z0) != 0:
        raise ArithmeticError("pairing condition misses the marked section")

    orthogonal = residual + ring.ideal([pairing_norm])
    orthogonal = ring.ideal(orthogonal.groebner_basis())
    orthogonal_degree = orthogonal.vector_space_dimension()
    print(f"orthogonal_degree={orthogonal_degree}", flush=True)

    nonorthogonal = residual.quotient(ring.ideal([pairing_norm]))
    nonorthogonal = ring.ideal(nonorthogonal.groebner_basis())
    nonorthogonal_degree = nonorthogonal.vector_space_dimension()
    print(f"nonorthogonal_degree={nonorthogonal_degree}", flush=True)

    section = ring.ideal([u - u0, v - v0, z - z0])
    orthogonal_residual = orthogonal.quotient(section)
    orthogonal_residual = ring.ideal(orthogonal_residual.groebner_basis())
    orthogonal_residual_degree = orthogonal_residual.vector_space_dimension()
    print(f"orthogonal_without_marked_degree={orthogonal_residual_degree}", flush=True)

    if (orthogonal_degree, orthogonal_residual_degree, nonorthogonal_degree) != (
        13, 12, 27
    ):
        raise ArithmeticError("Weil-pairing split did not have degrees 13,12,27")

    if check_intersections:
        reconstructs = orthogonal.intersection(nonorthogonal) == residual
        comaximal = ring.one() in orthogonal + nonorthogonal
        marked_reconstructs = (
            section.intersection(orthogonal_residual) == orthogonal
        )
        marked_comaximal = ring.one() in section + orthogonal_residual
        print(f"orbit_intersection_reconstructs={reconstructs}", flush=True)
        print(f"orbit_ideals_comaximal={comaximal}", flush=True)
        print(
            f"marked_intersection_reconstructs={marked_reconstructs}",
            flush=True,
        )
        print(f"marked_ideals_comaximal={marked_comaximal}", flush=True)
        if not all((reconstructs, comaximal, marked_reconstructs, marked_comaximal)):
            raise ArithmeticError("component reconstruction check failed")
    print(f"END base={label}", flush=True)


def main():
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--global", dest="global_base", action="store_true")
    parser.add_argument("--t", default="2")
    parser.add_argument("--skip-intersections", action="store_true")
    args = parser.parse_args()
    analyze(
        args.global_base,
        QQ(args.t),
        not args.skip_intersections,
    )


if __name__ == "__main__":
    main()
