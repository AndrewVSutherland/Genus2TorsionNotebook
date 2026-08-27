#!/usr/bin/env sage
"""Sage/Singular elimination for the unsigned contact-9 order-3 q-cover."""

import argparse

from sage.all import PolynomialRing, QQ


def system():
    function_field = QQ["t"].fraction_field()
    t = function_field.gen()
    x_ring = PolynomialRing(function_field, "x")
    x = x_ring.gen()

    r = 1 - t**2
    h_at_r = 1 - QQ(9)/2*r + QQ(63)/8*r**2 - QQ(105)/16*r**3
    contact_a = (t**9 - h_at_r)/r**4
    marked_h = (
        1 - QQ(9)/2*x + QQ(63)/8*x**2 - QQ(105)/16*x**3
        + contact_a*x**4
    )
    f, remainder = (marked_h**2 + (x - 1)**9).quo_rem(x**4)
    assert remainder == 0 and f.degree() == 5 and f[5] == 1

    ring = PolynomialRing(
        function_field, ("u", "v", "z"), order="degrevlex"
    )
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

    polynomial_over_ring = PolynomialRing(ring, "X")
    X = polynomial_over_ring.gen()
    q = X**2 + u*X + v
    f_over_ring = sum(coefficients[index]*X**index for index in range(6))
    resultant = q.resultant(f_over_ring)
    return function_field, ring, f, (a, b, c), equations, resultant


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--skip-resultant", action="store_true",
        help="saturate only by z and disc(q), for timing comparison",
    )
    parser.add_argument(
        "--method", choices=("quotient", "sequential", "product"),
        default="quotient",
    )
    parser.add_argument(
        "--skip-fglm", action="store_true",
        help="stop after the saturated degree computation",
    )
    parser.add_argument(
        "--skip-postfilters", action="store_true",
        help="after boundary quotients, do not saturate finite disc/resultant loci",
    )
    args = parser.parse_args()

    function_field, ring, f, triangular, equations, resultant = system()
    u, v, z = ring.gens()
    discriminant = u**2 - 4*v
    full_boundary = z*discriminant*resultant
    boundary = z*discriminant if args.skip_resultant else full_boundary
    raw = ring.ideal(equations)

    print("contact9 [3,18] unsigned q-cover Sage/Singular")
    print(f"coefficient_field={function_field}")
    print(f"f={f}")
    print(f"triangular_a={triangular[0]}")
    print(f"triangular_b={triangular[1]}")
    print(f"triangular_c={triangular[2]}")
    for index, equation in enumerate(equations, start=1):
        print(
            f"equation_{index}_degree={equation.degree()} "
            f"terms={len(equation.monomials())}"
        )
    print(f"raw_dimension={raw.dimension()}")
    print(f"boundary_degree={boundary.degree()}")

    if args.method == "quotient":
        boundary_component = ring.ideal([z, discriminant])
        saturated = raw
        for iteration in range(1, 5):
            saturated = saturated.quotient(boundary_component)
            print(
                f"boundary_quotient_iteration={iteration} "
                f"dimension={saturated.dimension()}"
            )
            if saturated.dimension() == 0:
                break
        if saturated.dimension() == 0 and not args.skip_postfilters:
            saturated, disc_power = saturated.saturation(
                ring.ideal([discriminant])
            )
            print(f"postquotient_discriminant_saturation_power={disc_power}")
            if not args.skip_resultant:
                saturated, resultant_power = saturated.saturation(
                    ring.ideal([resultant])
                )
                print(f"postquotient_resultant_saturation_power={resultant_power}")
    elif args.method == "sequential":
        saturated, z_power = raw.saturation(ring.ideal([z]))
        print(f"z_saturation_power={z_power}")
        saturated, disc_power = saturated.saturation(
            ring.ideal([discriminant])
        )
        print(f"discriminant_saturation_power={disc_power}")
        if not args.skip_resultant:
            saturated, resultant_power = saturated.saturation(
                ring.ideal([resultant])
            )
            print(f"resultant_saturation_power={resultant_power}")
    else:
        saturated, saturation_power = raw.saturation(ring.ideal([boundary]))
        print(f"saturation_power={saturation_power}")
    print(f"saturated_dimension={saturated.dimension()}")
    grevlex_basis = saturated.groebner_basis()
    saturated = ring.ideal(grevlex_basis)
    print(f"grevlex_basis_count={len(grevlex_basis)}")
    print(f"quotient_degree={saturated.vector_space_dimension()}")

    if args.skip_fglm:
        return
    lex_ring = PolynomialRing(
        function_field, ("u", "v", "z"), order="lex"
    )
    lex_basis = saturated.transformed_basis("fglm", lex_ring)
    print(f"lex_basis_count={len(lex_basis)}")
    for index, polynomial in enumerate(lex_basis, start=1):
        print(
            f"lex_basis_{index}_degree={polynomial.degree()} "
            f"terms={len(polynomial.monomials())}"
        )
    eliminants = [
        polynomial for polynomial in lex_basis
        if polynomial.degree(lex_ring.gen(0)) == 0
        and polynomial.degree(lex_ring.gen(1)) == 0
    ]
    print(f"z_eliminant_count={len(eliminants)}")
    for eliminant in eliminants:
        factorization = eliminant.factor()
        print(f"z_eliminant_degree={eliminant.degree(lex_ring.gen(2))}")
        print(
            "z_factor_degrees="
            + repr([
                (factor.degree(lex_ring.gen(2)), multiplicity)
                for factor, multiplicity in factorization
            ])
        )
        print(f"z_eliminant={eliminant}")
        for factor, multiplicity in factorization:
            print(
                f"z_factor degree={factor.degree(lex_ring.gen(2))} "
                f"multiplicity={multiplicity} polynomial={factor}"
            )


if __name__ == "__main__":
    main()
