#!/usr/bin/env sage
"""Multiplication-matrix factorization of the degree-40 contact-9 q-cover.

This avoids a lexicographic Groebner basis.  After removing the explicit
boundary component B=<z,u^2-4v> by two ideal quotients, it enumerates the
40 standard monomials of the grevlex quotient and builds regular
representation matrices by normal-form reduction.
"""

import argparse
from collections import deque

from sage.all import (
    GF,
    PolynomialRing,
    QQ,
    gcd,
    lcm,
    matrix,
    prod,
    vector,
)


LINEAR_FORMS = (
    (2, 3, 1),
    (1, 1, 1),
    (1, 2, 4),
    (3, 5, 7),
)


def qcover_system():
    parameter_ring = PolynomialRing(QQ, "t")
    function_field = parameter_ring.fraction_field()
    t = function_field(parameter_ring.gen())
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
    return parameter_ring, function_field, ring, f, equations


def open_quotient():
    parameter_ring, function_field, ring, f, equations = qcover_system()
    u, v, z = ring.gens()
    raw = ring.ideal(equations)
    boundary_component = ring.ideal([z, u**2 - 4*v])
    residual = raw.quotient(boundary_component)
    residual = residual.quotient(boundary_component)
    if residual.dimension() != 0:
        raise ArithmeticError("two boundary quotients did not give dimension zero")
    basis = tuple(residual.groebner_basis())
    residual = ring.ideal(basis)
    degree = residual.vector_space_dimension()
    if degree != 40:
        raise ArithmeticError(f"expected quotient degree 40, got {degree}")
    return parameter_ring, function_field, ring, f, residual, basis


def exponent_tuple(monomial):
    exponents = monomial.exponents()
    if len(exponents) != 1:
        raise ArithmeticError("expected a monomial")
    return tuple(exponents[0])


def divides(left, right):
    return all(a <= b for a, b in zip(left, right))


def standard_monomials(ring, groebner_basis):
    leading = tuple(
        exponent_tuple(polynomial.lm()) for polynomial in groebner_basis
    )
    zero = (0,)*ring.ngens()
    exponents = {zero}
    queue = deque([zero])
    while queue:
        current = queue.popleft()
        for index in range(ring.ngens()):
            candidate = list(current)
            candidate[index] += 1
            candidate = tuple(candidate)
            if candidate in exponents:
                continue
            if any(divides(lead, candidate) for lead in leading):
                continue
            exponents.add(candidate)
            queue.append(candidate)
            if len(exponents) > 1000:
                raise ArithmeticError("standard-monomial enumeration did not terminate")
    ordered_exponents = sorted(exponents, key=lambda item: (sum(item), item))
    monomials = tuple(
        prod(generator**power for generator, power in zip(ring.gens(), item))
        for item in ordered_exponents
    )
    return ordered_exponents, monomials, leading


def multiplication_matrix(
    function_field, groebner_basis, exponents, monomials, element
):
    index = {exponent: position for position, exponent in enumerate(exponents)}
    size = len(monomials)
    result = matrix(function_field, size, size)
    for column, monomial in enumerate(monomials):
        remainder = (element*monomial).reduce(groebner_basis)
        for exponent, coefficient in remainder.dict().items():
            exponent = tuple(exponent)
            if exponent not in index:
                raise ArithmeticError(
                    f"normal form contains nonstandard monomial {exponent}"
                )
            result[index[exponent], column] = coefficient
    return result


def evaluate_polynomial_matrix(function_field, polynomial, coordinate_matrices):
    size = coordinate_matrices[0].nrows()
    identity = matrix.identity(function_field, size)
    maxima = [0, 0, 0]
    for exponent in polynomial.dict():
        for index, power in enumerate(exponent):
            maxima[index] = max(maxima[index], power)
    powers = []
    for coordinate, maximum in zip(coordinate_matrices, maxima):
        coordinate_powers = [identity]
        for _index in range(maximum):
            coordinate_powers.append(coordinate_powers[-1]*coordinate)
        powers.append(coordinate_powers)
    result = matrix(function_field, size, size)
    for exponent, coefficient in polynomial.dict().items():
        term = identity
        for index, power in enumerate(exponent):
            term = term*powers[index][power]
        result += coefficient*term
    return result


def krylov_matrix(function_field, multiplication, identity_index):
    size = multiplication.nrows()
    current = vector(function_field, [int(i == identity_index) for i in range(size)])
    columns = []
    for _power in range(size):
        columns.append(current)
        current = multiplication*current
    return matrix(
        function_field, size, size,
        lambda row, column: columns[column][row],
    )


def berlekamp_massey(sequence):
    """Return recurrence C with s_n+sum(C_i*s_(n-i))=0."""
    field = sequence[0].parent()
    connection = [field(1)]
    previous = [field(1)]
    length = 0
    shift = 1
    previous_discrepancy = field(1)
    for position in range(len(sequence)):
        discrepancy = sequence[position]
        for index in range(1, length + 1):
            if index < len(connection):
                discrepancy += connection[index]*sequence[position - index]
        if discrepancy == 0:
            shift += 1
            continue
        old_connection = list(connection)
        scale = -discrepancy/previous_discrepancy
        required = len(previous) + shift
        if len(connection) < required:
            connection.extend([field(0)]*(required - len(connection)))
        for index, coefficient in enumerate(previous):
            connection[index + shift] += scale*coefficient
        if 2*length <= position:
            length = position + 1 - length
            previous = old_connection
            previous_discrepancy = discrepancy
            shift = 1
        else:
            shift += 1
    if len(connection) < length + 1:
        connection.extend([field(0)]*(length + 1 - len(connection)))
    return tuple(connection[:length + 1])


def wiedemann_minimal_polynomial(
    function_field, multiplication, identity_index, polynomial_ring
):
    size = multiplication.nrows()
    identity = vector(
        function_field, [int(index == identity_index) for index in range(size)]
    )
    functionals = (
        vector(function_field, [1]*size),
        vector(function_field, [index + 1 for index in range(size)]),
        vector(function_field, [(index + 1)**2 for index in range(size)]),
    )
    for functional_index, functional in enumerate(functionals):
        current = identity
        powers = []
        sequence = []
        for _index in range(2*size + 1):
            powers.append(current)
            sequence.append(functional.dot_product(current))
            current = multiplication*current
        recurrence = berlekamp_massey(sequence)
        degree = len(recurrence) - 1
        print(
            f"wiedemann_functional={functional_index} degree={degree}",
            flush=True,
        )
        if degree != size:
            continue
        relation = powers[degree]
        for index in range(1, degree + 1):
            relation += recurrence[index]*powers[degree - index]
        if relation:
            raise ArithmeticError("degree-40 scalar recurrence failed vector check")
        variable = polynomial_ring.gen()
        polynomial = variable**degree
        for index in range(1, degree + 1):
            polynomial += recurrence[index]*variable**(degree - index)
        return polynomial_ring(polynomial), powers[:size]
    return None, None


def marked_section(function_field):
    t = function_field.gen()
    denominator4 = (t + 1)**4
    u0 = (
        3*t**5 + QQ(15)/4*t**4 - 3*t**3 - QQ(147)/16*t**2
        - QQ(27)/4*t - QQ(27)/16
    )/denominator4
    v0 = (
        -2*t**5 - QQ(47)/16*t**4 + QQ(1)/4*t**3
        + QQ(7)/2*t**2 + QQ(11)/4*t + QQ(11)/16
    )/denominator4
    z0 = (
        4*t**10 + 14*t**9 + QQ(81)/4*t**8 + QQ(27)/2*t**7
        + QQ(9)/8*t**6 - QQ(9)/2*t**5 - QQ(183)/64*t**4
        - QQ(3)/8*t**3 + QQ(9)/32*t**2 + QQ(1)/8*t + QQ(1)/64
    )/(t + 1)**8
    return u0, v0, z0


def weil_pairing_norm(function_field, ring, f, residual):
    t = function_field.gen()
    u, v, z = ring.gens()
    coefficients = [ring(f[index]) for index in range(6)]
    a = (3*u + z)/2
    b = (z*coefficients[4] + 3*u**2 + 3*v - a**2)/2
    c = (z*coefficients[3] + u**3 + 6*u*v - 2*a*b)/2

    u0, v0, z0 = marked_section(function_field)
    s0 = (
        2*(t + QQ(1)/2)**2*(t**3 + QQ(3)/4*t**2 - QQ(1)/4)
        / (t + 1)**4
    )
    if s0**2 != z0:
        raise ArithmeticError("marked square root does not square to z0")
    a0 = (3*u0 + z0)/2
    b0 = (
        z0*function_field(f[4]) + 3*u0**2 + 3*v0 - a0**2
    )/2
    c0 = (
        z0*function_field(f[3]) + u0**3 + 6*u0*v0 - 2*a0*b0
    )/2

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
    value_on_candidate = q.resultant(difference)
    pairing_equation = lam**2*value_on_marked - value_on_candidate
    relation = lam**2 - ring(z/z0)
    _quotient, remainder = pairing_equation.quo_rem(relation)
    if remainder.degree() > 1:
        raise ArithmeticError("pairing equation did not reduce linearly in lambda")
    constant = ring(remainder[0])
    linear = ring(remainder[1])
    norm = constant**2 - ring(z/z0)*linear**2
    return residual.reduce(norm), (constant, linear)


def coordinate_polynomial(polynomial_ring, krylov, target):
    coefficients = krylov.solve_right(target)
    variable = polynomial_ring.gen()
    return polynomial_ring(
        sum(coefficient*variable**index
            for index, coefficient in enumerate(coefficients))
    )


def clear_parameter_denominators(parameter_ring, polynomial):
    denominators = [
        parameter_ring(coefficient.denominator())
        for coefficient in polynomial.list()
    ]
    common = parameter_ring(1)
    for denominator in denominators:
        common = lcm(common, denominator)
    bivariate = PolynomialRing(QQ, ("t", "L"))
    tvar, Lvar = bivariate.gens()
    cleared = bivariate(0)
    for index, coefficient in enumerate(polynomial.list()):
        numerator = parameter_ring(coefficient.numerator())
        denominator = parameter_ring(coefficient.denominator())
        multiplier = common // denominator
        coefficient_polynomial = numerator*multiplier
        cleared += sum(
            QQ(coefficient_polynomial[degree])*tvar**degree
            for degree in range(coefficient_polynomial.degree() + 1)
        )*Lvar**index
    rational_content = gcd(
        [coefficient for coefficient in cleared.coefficients()]
    )
    if rational_content:
        cleared /= rational_content
    return cleared


def reduce_function_at_4_mod_7(value, parameter):
    field = GF(7)

    def evaluate_polynomial(poly):
        result = field(0)
        for index in range(poly.degree() + 1):
            result += field(poly[index])*field(4)**index
        return result

    numerator = parameter(value.numerator())
    denominator = parameter(value.denominator())
    denominator_value = evaluate_polynomial(denominator)
    if denominator_value == 0:
        return None
    return evaluate_polynomial(numerator)/denominator_value


def reduce_polynomial_at_4_mod_7(polynomial, parameter):
    finite_ring = PolynomialRing(GF(7), "L")
    coefficients = []
    for coefficient in polynomial.list():
        reduced = reduce_function_at_4_mod_7(coefficient, parameter)
        if reduced is None:
            return None
        coefficients.append(reduced)
    return finite_ring(coefficients)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--forms", default="all",
        help="'all' or a comma-separated list of indices starting at 0",
    )
    parser.add_argument("--skip-coordinates", action="store_true")
    parser.add_argument(
        "--algorithm", choices=("wiedemann", "charpoly"), default="wiedemann"
    )
    parser.add_argument(
        "--export-magma",
        help="write the first requested multiplication matrix as Magma source and stop",
    )
    parser.add_argument("--weil-rank", action="store_true")
    args = parser.parse_args()

    print("stage=build_open_quotient", flush=True)
    (
        parameter_ring, function_field, ring, f, residual, groebner_basis
    ) = open_quotient()
    print("stage=open_quotient_complete", flush=True)
    u, v, z = ring.gens()
    exponents, monomials, leading = standard_monomials(ring, groebner_basis)
    print("stage=standard_monomials_complete", flush=True)
    print("contact9 [3,18] q-cover multiplication matrices")
    print(f"quotient_degree={residual.vector_space_dimension()}")
    print(f"groebner_basis_count={len(groebner_basis)}")
    print(f"leading_exponents={leading}")
    print(f"standard_monomial_count={len(monomials)}")
    print(f"standard_exponents={exponents}")
    if len(monomials) != 40:
        raise ArithmeticError("standard monomial count is not 40")
    identity_index = exponents.index((0, 0, 0))

    if args.weil_rank:
        print("stage=weil_pairing_norm_start", flush=True)
        pairing_norm, pairing_remainder = weil_pairing_norm(
            function_field, ring, f, residual
        )
        print(
            f"pairing_norm_degree={pairing_norm.degree()} "
            f"pairing_norm_terms={len(pairing_norm.monomials())}",
            flush=True,
        )
        print(
            f"pairing_remainder_degrees="
            f"({pairing_remainder[0].degree()},{pairing_remainder[1].degree()})",
            flush=True,
        )
        print("stage=weil_coordinate_matrices_start", flush=True)
        coordinate_matrices = []
        for name, generator in zip(("u", "v", "z"), ring.gens()):
            coordinate_matrices.append(
                multiplication_matrix(
                    function_field, groebner_basis, exponents, monomials,
                    generator,
                )
            )
            print(f"stage=weil_coordinate_matrix_{name}_complete", flush=True)
        print("stage=weil_matrix_evaluation_start", flush=True)
        weil_matrix = evaluate_polynomial_matrix(
            function_field, pairing_norm, coordinate_matrices
        )
        print("stage=weil_matrix_complete", flush=True)
        print("stage=weil_rank_start", flush=True)
        weil_rank = weil_matrix.rank()
        print(f"weil_rank={weil_rank}", flush=True)
        print(f"weil_nullity={40 - weil_rank}", flush=True)
        return

    if args.forms == "all":
        form_indices = range(len(LINEAR_FORMS))
    else:
        form_indices = [int(index) for index in args.forms.split(",")]

    coordinate_matrices = None
    selected = None
    for form_index in form_indices:
        coefficients = LINEAR_FORMS[form_index]
        linear_form = sum(
            coefficient*generator
            for coefficient, generator in zip(coefficients, ring.gens())
        )
        print(f"stage=form_{form_index}_matrix_start", flush=True)
        multiplication = multiplication_matrix(
            function_field, groebner_basis, exponents, monomials, linear_form
        )
        print(f"stage=form_{form_index}_matrix_complete", flush=True)
        if args.export_magma:
            with open(args.export_magma, "w", encoding="ascii") as output:
                output.write("Q := Rationals();\n")
                output.write("K<t> := FunctionField(Q);\n")
                output.write("PL<L> := PolynomialRing(K);\n")
                output.write("M := Matrix(K, 40, 40, [\n")
                entries = [
                    repr(multiplication[row, column])
                    for row in range(40) for column in range(40)
                ]
                for start in range(0, len(entries), 4):
                    suffix = "," if start + 4 < len(entries) else ""
                    output.write("  " + ", ".join(entries[start:start + 4]) + suffix + "\n")
                output.write("]);\n")
                output.write(
                    "marked_u := (3*t^5+15/4*t^4-3*t^3-147/16*t^2"
                    "-27/4*t-27/16)/(t+1)^4;\n"
                )
                output.write(
                    "marked_v := (-2*t^5-47/16*t^4+1/4*t^3+7/2*t^2"
                    "+11/4*t+11/16)/(t+1)^4;\n"
                )
                output.write(
                    "marked_z := (4*t^10+14*t^9+81/4*t^8+27/2*t^7"
                    "+9/8*t^6-9/2*t^5-183/64*t^4-3/8*t^3"
                    "+9/32*t^2+1/8*t+1/64)/(t+1)^8;\n"
                )
                output.write(
                    f"marked_value := {coefficients[0]}*marked_u"
                    f"+{coefficients[1]}*marked_v"
                    f"+{coefficients[2]}*marked_z;\n"
                )
            print(f"exported_magma_matrix={args.export_magma}", flush=True)
            return
        characteristic_ring = PolynomialRing(function_field, "L")
        if args.algorithm == "charpoly":
            print(f"stage=form_{form_index}_charpoly_start", flush=True)
            characteristic = multiplication.charpoly("L")
            print(f"stage=form_{form_index}_charpoly_complete", flush=True)
        else:
            print(f"stage=form_{form_index}_wiedemann_start", flush=True)
            characteristic, _power_vectors = wiedemann_minimal_polynomial(
                function_field, multiplication, identity_index,
                characteristic_ring,
            )
            print(f"stage=form_{form_index}_wiedemann_complete", flush=True)
            if characteristic is None:
                continue
        characteristic_ring = characteristic.parent()
        derivative_gcd = characteristic.gcd(characteristic.derivative())
        squarefree = derivative_gcd.degree() == 0
        print(
            f"form_index={form_index} coefficients={coefficients} "
            f"characteristic_squarefree={squarefree}"
        )
        if not squarefree:
            continue
        marked_u, marked_v, marked_z = marked_section(function_field)
        marked_value = (
            coefficients[0]*marked_u + coefficients[1]*marked_v
            + coefficients[2]*marked_z
        )
        marked_linear = characteristic_ring.gen() - marked_value
        marked_quotient, marked_remainder = characteristic.quo_rem(
            marked_linear
        )
        print(f"marked_linear_root={marked_value}")
        print(f"marked_linear_remainder={marked_remainder}")
        if marked_remainder:
            raise ArithmeticError("known marked section is not a polynomial root")
        print(f"stage=form_{form_index}_factor_start", flush=True)
        residual_factorization = marked_quotient.factor()
        factorization = [(marked_linear, 1)] + list(residual_factorization)
        print(f"stage=form_{form_index}_factor_complete", flush=True)
        print(f"selected_form_index={form_index}")
        print(f"characteristic_degree={characteristic.degree()}")
        print(f"characteristic_squarefree={derivative_gcd.degree() == 0}")
        print(
            "factor_degrees="
            + repr([
                (factor.degree(), multiplicity)
                for factor, multiplicity in factorization
            ])
        )
        print(f"characteristic_polynomial={characteristic}")
        for factor_index, (factor, multiplicity) in enumerate(factorization):
            cleared = clear_parameter_denominators(parameter_ring, factor)
            print(
                f"factor_index={factor_index} degree={factor.degree()} "
                f"multiplicity={multiplicity} "
                f"bidegree_t_L=({cleared.degree(cleared.parent().gen(0))},"
                f"{cleared.degree(cleared.parent().gen(1))}) "
                f"terms={len(cleared.monomials())}"
            )
            print(f"factor_{factor_index}={factor}")
            print(f"factor_{factor_index}_cleared={cleared}")
            reduction = reduce_polynomial_at_4_mod_7(
                factor, parameter_ring
            )
            print(f"factor_{factor_index}_t4_mod7={reduction}")
            if reduction is not None:
                print(
                    f" factor_at_branch0_L1={reduction(GF(7)(1))} "
                    f"factor_at_branch1_L2={reduction(GF(7)(2))}"
                )
        selected = (
            coefficients, multiplication, characteristic,
            characteristic_ring, factorization,
        )
        break

    if selected is None:
        print("NO_SEPARATING_LINEAR_FORM_IN_DECLARED_LIST")
        return
    if args.skip_coordinates:
        return

    if coordinate_matrices is None:
        coordinate_matrices = [
            multiplication_matrix(
                function_field, groebner_basis, exponents, monomials, generator
            )
            for generator in (u, v, z)
        ]
    coefficients, multiplication, characteristic, poly_ring, _factors = selected
    print("stage=coordinate_krylov_start", flush=True)
    krylov = krylov_matrix(function_field, multiplication, identity_index)
    krylov_rank = krylov.rank()
    print(f"coordinate_krylov_rank={krylov_rank}", flush=True)
    if krylov_rank != 40:
        raise ArithmeticError("squarefree characteristic polynomial but Krylov rank < 40")
    identity_vector = vector(
        function_field,
        [int(index == identity_index) for index in range(40)],
    )
    coordinate_polynomials = []
    for name, coordinate_matrix in zip(("u", "v", "z"), coordinate_matrices):
        polynomial = coordinate_polynomial(
            poly_ring, krylov, coordinate_matrix*identity_vector
        )
        coordinate_polynomials.append(polynomial)
        print(f"coordinate_{name}={polynomial}")

    relation = (
        coefficients[0]*coordinate_polynomials[0]
        + coefficients[1]*coordinate_polynomials[1]
        + coefficients[2]*coordinate_polynomials[2]
        - poly_ring.gen()
    ) % characteristic
    print(f"linear_form_coordinate_relation={relation}")

    for name, polynomial in zip(("u", "v", "z"), coordinate_polynomials):
        reduction = reduce_polynomial_at_4_mod_7(polynomial, parameter_ring)
        print(f"coordinate_{name}_t4_mod7={reduction}")
        if reduction is not None:
            print(
                f" coordinate_{name}_branch0_L1={reduction(GF(7)(1))} "
                f"coordinate_{name}_branch1_L2={reduction(GF(7)(2))}"
            )


if __name__ == "__main__":
    main()
