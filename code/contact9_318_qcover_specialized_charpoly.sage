#!/usr/bin/env sage
"""Exact rational-fiber diagnostics for the contact-9 order-3 q-cover.

This is deliberately separate from the global Q(t) calculation.  For each
specified rational value of t it constructs the same double ideal quotient,
builds multiplication by ell=2*u+3*v+z in the degree-40 algebra, and factors
the characteristic polynomial over Q.
"""

import argparse

from sage.all import PolynomialRing, QQ, matrix


def build_fiber(parameter):
    x_ring = PolynomialRing(QQ, "x")
    x = x_ring.gen()
    t = QQ(parameter)

    r = 1 - t**2
    if r == 0:
        raise ValueError("t=+-1 is excluded from this presentation")
    h_at_r = 1 - QQ(9)/2*r + QQ(63)/8*r**2 - QQ(105)/16*r**3
    contact_a = (t**9 - h_at_r)/r**4
    marked_h = (
        1 - QQ(9)/2*x + QQ(63)/8*x**2 - QQ(105)/16*x**3
        + contact_a*x**4
    )
    f, remainder = (marked_h**2 + (x - 1)**9).quo_rem(x**4)
    if remainder != 0 or f.degree() != 5 or f[5] != 1:
        raise ArithmeticError("bad contact-9 quintic")

    ring = PolynomialRing(QQ, ("u", "v", "z"), order="degrevlex")
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
    return ring, f, equations


def marked_coordinates(t):
    d4 = (t + 1)**4
    u0 = (
        3*t**5 + QQ(15)/4*t**4 - 3*t**3 - QQ(147)/16*t**2
        - QQ(27)/4*t - QQ(27)/16
    )/d4
    v0 = (
        -2*t**5 - QQ(47)/16*t**4 + QQ(1)/4*t**3
        + QQ(7)/2*t**2 + QQ(11)/4*t + QQ(11)/16
    )/d4
    z0 = (
        4*t**10 + 14*t**9 + QQ(81)/4*t**8 + QQ(27)/2*t**7
        + QQ(9)/8*t**6 - QQ(9)/2*t**5 - QQ(183)/64*t**4
        - QQ(3)/8*t**3 + QQ(9)/32*t**2 + QQ(1)/8*t
        + QQ(1)/64
    )/(t + 1)**8
    return u0, v0, z0


def coefficient_vector(polynomial, monomial_basis):
    return [
        polynomial.monomial_coefficient(monomial)
        for monomial in monomial_basis
    ]


def analyze(parameter):
    t = QQ(parameter)
    ring, f, equations = build_fiber(t)
    u, v, z = ring.gens()
    raw = ring.ideal(equations)
    boundary = ring.ideal([z, u**2 - 4*v])

    print(f"BEGIN t={t}", flush=True)
    print(f"curve_squarefree={f.gcd(f.derivative()).degree() == 0}", flush=True)
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

    basis = residual.groebner_basis()
    residual = ring.ideal(basis)
    degree = residual.vector_space_dimension()
    print(f"groebner_basis_count={len(basis)}", flush=True)
    print(f"quotient_degree={degree}", flush=True)
    if degree != 40:
        raise ArithmeticError(f"expected degree 40, got {degree}")

    monomials = residual.normal_basis()
    ell = 2*u + 3*v + z
    columns = []
    for monomial in monomials:
        normal_form = residual.reduce(ell*monomial)
        coordinates = coefficient_vector(normal_form, monomials)
        reconstructed = sum(
            coordinates[index]*monomials[index]
            for index in range(degree)
        )
        if reconstructed != normal_form:
            raise ArithmeticError("normal-form coordinate extraction failed")
        columns.append(coordinates)
    multiplication = matrix(
        QQ, degree, degree,
        lambda row, column: columns[column][row],
    )
    characteristic = multiplication.charpoly("L")
    squarefree = characteristic.gcd(characteristic.derivative()).degree() == 0
    print(f"charpoly_degree={characteristic.degree()}", flush=True)
    print(f"charpoly_squarefree={squarefree}", flush=True)

    one = matrix(QQ, degree, 1, coefficient_vector(ring.one(), monomials))
    krylov_columns = []
    current = one
    for _ in range(degree):
        krylov_columns.append(list(current.column(0)))
        current = multiplication*current
    krylov = matrix(
        QQ, degree, degree,
        lambda row, column: krylov_columns[column][row],
    )
    print(f"ell_krylov_rank={krylov.rank()}", flush=True)

    factorization = characteristic.factor()
    pattern = [
        (factor.degree(), multiplicity)
        for factor, multiplicity in factorization
    ]
    print(f"charpoly_factor_pattern={pattern}", flush=True)
    u0, v0, z0 = marked_coordinates(t)
    ell0 = 2*u0 + 3*v0 + z0
    print(f"marked_coordinates={(u0, v0, z0)}", flush=True)
    print(f"marked_ell={ell0}", flush=True)
    print(f"marked_charpoly_value={characteristic(ell0)}", flush=True)
    for index, (factor, multiplicity) in enumerate(factorization):
        print(
            f"factor_{index}_degree={factor.degree()} "
            f"multiplicity={multiplicity} marked_value={factor(ell0)}",
            flush=True,
        )
    print(f"END t={t}", flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--values", default="2,3,4",
        help="comma-separated rational values, for example 2,3/2,4",
    )
    args = parser.parse_args()
    for value in args.values.split(","):
        analyze(QQ(value.strip()))


if __name__ == "__main__":
    main()
