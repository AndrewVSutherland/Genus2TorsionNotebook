#!/usr/bin/env sage
"""Modular multiplication-polynomial diagnostics for the [3,18] q-cover.

For one odd prime p, construct the unsigned q-cover over F_p(t), remove the
boundary by the same double ideal quotient used over Q(t), and study the
degree-40 residual algebra with the multiplication operator

    ell = z + 2*u + 3*v.

Run each prime under an external hard timeout.
"""

import argparse

from sage.all import GF, PolynomialRing, matrix


def build_system(prime):
    constants = GF(prime)
    parameter_ring = PolynomialRing(constants, "t")
    function_field = parameter_ring.fraction_field()
    t = function_field.gen()
    x_ring = PolynomialRing(function_field, "x")
    x = x_ring.gen()

    def c(value):
        return function_field(value)

    r = 1 - t**2
    h_at_r = 1 - c(9)/2*r + c(63)/8*r**2 - c(105)/16*r**3
    contact_a = (t**9 - h_at_r)/r**4
    marked_h = (
        1 - c(9)/2*x + c(63)/8*x**2 - c(105)/16*x**3
        + contact_a*x**4
    )
    f, remainder = (marked_h**2 + (x - 1)**9).quo_rem(x**4)
    if remainder != 0 or f.degree() != 5 or f[5] != 1:
        raise ValueError("bad reduction of the contact-9 quintic")

    ring = PolynomialRing(
        function_field, ("u", "v", "z"), order="degrevlex"
    )
    u, v, z = ring.gens()
    coefficients = [ring(f[index]) for index in range(6)]
    a = (3*u + z)/2
    b = (z*coefficients[4] + 3*u**2 + 3*v - a**2)/2
    cc = (z*coefficients[3] + u**3 + 6*u*v - 2*a*b)/2
    equations = (
        b**2 + 2*a*cc - z*coefficients[2] - 3*(u**2*v + v**2),
        2*b*cc - z*coefficients[1] - 3*u*v**2,
        cc**2 - z*coefficients[0] - v**3,
    )
    return function_field, ring, equations


def coefficient_vector(polynomial, monomial_basis):
    return [polynomial.monomial_coefficient(monomial)
            for monomial in monomial_basis]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--fglm", action="store_true")
    args = parser.parse_args()
    p = args.prime
    if p == 2:
        raise ValueError("characteristic 2 is excluded")

    function_field, ring, equations = build_system(p)
    u, v, z = ring.gens()
    raw = ring.ideal(equations)
    boundary_ideal = ring.ideal([z, u**2 - 4*v])

    print(f"QCover modular multiplication diagnostic p={p}", flush=True)
    print(f"coefficient_field={function_field}", flush=True)
    print(
        "equation_summaries="
        + repr([(g.degree(), len(g.monomials())) for g in equations]),
        flush=True,
    )
    print(f"raw_dimension={raw.dimension()}", flush=True)

    residual = raw
    for iteration in (1, 2):
        residual = residual.quotient(boundary_ideal)
        print(
            f"boundary_quotient_iteration={iteration} "
            f"dimension={residual.dimension()}",
            flush=True,
        )
    if residual.dimension() != 0:
        raise ValueError("double quotient is not zero-dimensional")

    grevlex_basis = residual.groebner_basis()
    residual = ring.ideal(grevlex_basis)
    degree = residual.vector_space_dimension()
    print(f"grevlex_basis_count={len(grevlex_basis)}", flush=True)
    print(f"quotient_degree={degree}", flush=True)
    if degree != 40:
        raise ValueError(f"bad modular residual degree {degree}, expected 40")

    monomial_basis = residual.normal_basis()
    if len(monomial_basis) != degree:
        raise ValueError("normal basis size disagrees with quotient degree")
    print(
        "normal_basis_degree_profile="
        + repr(sorted(m.degree() for m in monomial_basis)),
        flush=True,
    )

    ell = z + 2*u + 3*v
    columns = []
    for monomial in monomial_basis:
        normal_form = residual.reduce(ell*monomial)
        coordinates = coefficient_vector(normal_form, monomial_basis)
        reconstructed = sum(
            coordinates[index]*monomial_basis[index]
            for index in range(degree)
        )
        if reconstructed != normal_form:
            raise ValueError("normal-form coordinate extraction failed")
        columns.append(coordinates)
    multiplication = matrix(
        function_field, degree, degree,
        lambda row, column: columns[column][row],
    )
    characteristic = multiplication.charpoly("W")
    print(f"charpoly_degree={characteristic.degree()}", flush=True)
    squarefree = characteristic.gcd(characteristic.derivative()).degree() == 0
    print(f"charpoly_squarefree={squarefree}", flush=True)

    one_coordinates = coefficient_vector(ring.one(), monomial_basis)
    vector = matrix(function_field, degree, 1, one_coordinates)
    krylov_columns = []
    for _ in range(degree):
        krylov_columns.append(list(vector.column(0)))
        vector = multiplication*vector
    krylov = matrix(
        function_field, degree, degree,
        lambda row, column: krylov_columns[column][row],
    )
    cyclic_rank = krylov.rank()
    print(f"ell_krylov_rank={cyclic_rank}", flush=True)
    print(f"ell_separating={squarefree and cyclic_rank == degree}", flush=True)

    factors = characteristic.factor()
    factor_pattern = [
        (factor.degree(), multiplicity) for factor, multiplicity in factors
    ]
    print(f"charpoly_factor_pattern={factor_pattern}", flush=True)
    print(
        "charpoly_coefficient_t_degrees="
        + repr([
            (coefficient.numerator().degree(), coefficient.denominator().degree())
            if coefficient else (-1, -1)
            for coefficient in characteristic.list()
        ]),
        flush=True,
    )

    if args.fglm:
        lex_ring = PolynomialRing(
            function_field, ("u", "v", "z"), order="lex"
        )
        lex_basis = residual.transformed_basis("fglm", lex_ring)
        print(f"lex_basis_count={len(lex_basis)}", flush=True)
        print(
            "lex_basis_summaries="
            + repr([
                (g.degree(), len(g.monomials())) for g in lex_basis
            ]),
            flush=True,
        )
    print("DONE", flush=True)


if __name__ == "__main__":
    main()
