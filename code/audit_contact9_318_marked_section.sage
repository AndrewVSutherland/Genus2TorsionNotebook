#!/usr/bin/env sage
"""Extract and remove the tautological 3*D9 section of the [3,18] q-cover."""

import argparse
import sys

from sage.all import GF, HyperellipticCurve, PolynomialRing, QQ


def family():
    parameter_ring = PolynomialRing(QQ, "t")
    t_polynomial = parameter_ring.gen()
    field = parameter_ring.fraction_field()
    t = field(t_polynomial)
    x_ring = PolynomialRing(field, "x")
    x = x_ring.gen()
    root = 1 - t**2
    h_at_root = (
        1 - QQ(9) / 2 * root + QQ(63) / 8 * root**2
        - QQ(105) / 16 * root**3
    )
    contact_a = (t**9 - h_at_root) / root**4
    marked_h = (
        1 - QQ(9) / 2 * x + QQ(63) / 8 * x**2
        - QQ(105) / 16 * x**3 + contact_a * x**4
    )
    f, remainder = (marked_h**2 + (x - 1) ** 9).quo_rem(x**4)
    assert remainder == 0 and f.degree() == 5 and f[5] == 1
    return parameter_ring, field, t, x_ring, x, marked_h, f


def qcover_system(field, f):
    ring = PolynomialRing(field, ("u", "v", "z"), order="degrevlex")
    u, v, z = ring.gens()
    coefficients = [ring(f[index]) for index in range(6)]
    cubic_a = (3 * u + z) / 2
    cubic_b = (
        z * coefficients[4] + 3 * u**2 + 3 * v - cubic_a**2
    ) / 2
    cubic_c = (
        z * coefficients[3] + u**3 + 6 * u * v - 2 * cubic_a * cubic_b
    ) / 2
    equations = (
        cubic_b**2 + 2 * cubic_a * cubic_c
        - z * coefficients[2] - 3 * (u**2 * v + v**2),
        2 * cubic_b * cubic_c - z * coefficients[1] - 3 * u * v**2,
        cubic_c**2 - z * coefficients[0] - v**3,
    )
    polynomial_ring = PolynomialRing(ring, "X")
    X = polynomial_ring.gen()
    q = X**2 + u * X + v
    f_over_ring = sum(coefficients[index] * X**index for index in range(6))
    resultant = q.resultant(f_over_ring)
    return ring, (cubic_a, cubic_b, cubic_c), equations, resultant


def specialize(value, parameter, target):
    value = value.parent().fraction_field()(value)
    numerator = value.numerator()(parameter)
    denominator = value.denominator()(parameter)
    return target(numerator) / target(denominator)


def main():
    sys.stdout.reconfigure(line_buffering=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("--split", action="store_true")
    parser.add_argument("--section-only", action="store_true")
    parser.add_argument("--postfilters", action="store_true")
    parser.add_argument("--strong-check", action="store_true")
    args = parser.parse_args()

    parameter_ring, field, t, x_ring, x, marked_h, f = family()
    curve = HyperellipticCurve(f)
    jacobian = curve.jacobian()(field)
    zero = jacobian(0)
    marked = jacobian([x - 1, x_ring(marked_h(1))])
    marked_three = 3 * marked
    assert 9 * marked == zero and marked_three != zero

    q = x_ring(marked_three[0])
    mumford_v = x_ring(marked_three[1])
    u0, v0 = field(q[1]), field(q[0])

    ring, triangular, equations, resultant = qcover_system(field, f)
    u, v, z = ring.gens()
    z_ring = PolynomialRing(field, "Z")
    Z = z_ring.gen()
    section_equations = [
        z_ring(equation(u0, v0, Z)) for equation in equations
    ]
    common = section_equations[0].gcd(section_equations[1]).gcd(
        section_equations[2]
    ).monic()
    assert common.degree() == 1
    z0 = -common[0]

    cubic_a0 = field(triangular[0](u0, v0, z0))
    cubic_b0 = field(triangular[1](u0, v0, z0))
    cubic_c0 = field(triangular[2](u0, v0, z0))
    G = x**3 + cubic_a0 * x**2 + cubic_b0 * x + cubic_c0

    square_root_z = (
        2 * (t + QQ(1) / 2) ** 2
        * (t**3 + QQ(3) / 4 * t**2 - QQ(1) / 4)
        / (t + 1) ** 4
    )
    assert square_root_z**2 == z0
    contact_H = G / square_root_z
    relation_sign = None
    if contact_H % q == mumford_v:
        relation_sign = 1
    elif contact_H % q == -mumford_v:
        relation_sign = -1
    assert relation_sign is not None

    identity = G**2 - z0 * f - q**3
    residuals = tuple(equation(u0, v0, z0) for equation in equations)
    section_discriminant = u0**2 - 4 * v0
    section_resultant = resultant(u0, v0, z0)
    assert identity == 0 and not any(residuals)
    assert z0 != 0 and section_discriminant != 0 and section_resultant != 0

    print("CONTACT9_MARKED_SECTION")
    print(f"q={q}")
    print(f"u0={u0}")
    print(f"v0={v0}")
    print(f"z0={z0}")
    print(f"z0_factor={z0.factor()}")
    print(f"sqrt_z={square_root_z}")
    print(f"G={G}")
    print(f"contact_H=G/sqrt_z={contact_H}")
    print(f"mumford_v={mumford_v}")
    print(f"contact_H_mod_q_sign={relation_sign}")
    print(f"identity_residual={identity}")
    print(f"coefficient_residuals={residuals}")
    print(f"section_disc_q={section_discriminant}")
    print(f"section_resultant_q_f={section_resultant}")
    print("section_open=True")

    field7 = GF(7)
    t4_u = specialize(u0, QQ(4), field7)
    t4_v = specialize(v0, QQ(4), field7)
    t4_z = specialize(z0, QQ(4), field7)
    print(f"mod7_t4_q=x^2+{t4_u}*x+{t4_v}")
    print(f"mod7_t4_z={t4_z}")
    assert (t4_u, t4_v) == (field7(1), field7(4))
    if args.section_only:
        return

    raw = ring.ideal(equations)
    discriminant = u**2 - 4 * v
    boundary_component = ring.ideal([z, discriminant])
    saturated = raw
    quotient_iterations = 0
    while saturated.dimension() != 0 and quotient_iterations < 4:
        saturated = saturated.quotient(boundary_component)
        quotient_iterations += 1
        print(
            f"boundary_quotient_iteration={quotient_iterations} "
            f"dimension={saturated.dimension()}"
        )
    assert saturated.dimension() == 0
    disc_power = None
    resultant_power = None
    if args.postfilters:
        saturated, disc_power = saturated.saturation(ring.ideal([discriminant]))
        saturated, resultant_power = saturated.saturation(ring.ideal([resultant]))
    saturated = ring.ideal(saturated.groebner_basis())
    total_degree = saturated.vector_space_dimension()
    print(f"total_cover_degree={total_degree}")

    section = ring.ideal([u - u0, v - v0, z - z0])
    assert all(polynomial.reduce(section) == 0 for polynomial in saturated.gens())
    residual = saturated.quotient(section)
    residual = ring.ideal(residual.groebner_basis())
    residual_degree = residual.vector_space_dimension()
    print(f"residual_degree={residual_degree}")
    reconstruction_ok = None
    comaximal = None
    if args.strong_check:
        reconstructed = section.intersection(residual)
        reconstruction_ok = reconstructed == saturated
        comaximal = ring.one() in section + residual

    print(f"boundary_quotient_iterations={quotient_iterations}")
    print(f"discriminant_saturation_power={disc_power}")
    print(f"resultant_saturation_power={resultant_power}")
    print(f"section_degree={section.vector_space_dimension()}")
    print(f"intersection_reconstructs_cover={reconstruction_ok}")
    print(f"section_residual_comaximal={comaximal}")
    assert total_degree == 40 and residual_degree == 39
    if args.strong_check:
        assert reconstruction_ok and comaximal

    if args.split:
        lex_ring = PolynomialRing(field, ("u", "v", "z"), order="lex")
        lex_basis = residual.transformed_basis("fglm", lex_ring)
        eliminants = [
            polynomial for polynomial in lex_basis
            if polynomial.degree(lex_ring.gen(0)) == 0
            and polynomial.degree(lex_ring.gen(1)) == 0
        ]
        print(f"residual_lex_basis_count={len(lex_basis)}")
        print(f"residual_z_eliminant_count={len(eliminants)}")
        for eliminant in eliminants:
            factors = eliminant.factor()
            print(f"residual_z_eliminant_degree={eliminant.degree(lex_ring.gen(2))}")
            print(
                "residual_z_factor_degrees="
                + repr([
                    (factor.degree(lex_ring.gen(2)), multiplicity)
                    for factor, multiplicity in factors
                ])
            )


if __name__ == "__main__":
    main()
