#!/usr/bin/env sage
"""Resolve and search the order-3 cover above the contact-9 t=4 mod 7 chart.

For q=x^2+u*x+v and H=h3*x^3+h2*x^2+h1*x+h0, the identity

    H^2 - f = h3^2*q^3

gives a degree-two divisor D=(q,H mod q) with 3*D=0, provided the usual
open conditions hold.  This script exhausts the special fiber at t=4,
checks its tangent geometry, distinguishes the independent classes in the
generalized Jacobian of the nodal special curve, and performs a bounded
one-prime Hensel/rational-reconstruction search above those branches.
"""

import argparse
from collections import Counter
from itertools import product
from math import isqrt

from sage.all import (
    EllipticCurve,
    GF,
    HyperellipticCurve,
    Integer,
    PolynomialRing,
    QQ,
    ZZ,
    gcd,
    matrix,
    vector,
)


PRIME = 7
VARIABLE_NAMES = ("u", "v", "h0", "h1", "h2", "h3")


def contact9_polynomials(base_ring, t):
    """Return (r,a,h,f) for the rational-root contact-9 family."""
    polynomial_ring = PolynomialRing(base_ring, "x")
    x = polynomial_ring.gen()
    r = 1 - t**2
    h_at_r_without_a = (
        1 - base_ring(9)/2*r + base_ring(63)/8*r**2
        - base_ring(105)/16*r**3
    )
    a = (t**9 - h_at_r_without_a)/r**4
    h = (
        1 - base_ring(9)/2*x + base_ring(63)/8*x**2
        - base_ring(105)/16*x**3 + a*x**4
    )
    numerator = h**2 + (x - 1)**9
    f, remainder = numerator.quo_rem(x**4)
    if remainder:
        raise ArithmeticError("contact-9 division by x^4 failed")
    return r, a, h, f


def fiber_equations(f, values):
    """Six non-leading coefficients of H^2-f-h3^2*q^3."""
    u, v, h0, h1, h2, h3 = values
    polynomial_ring = f.parent()
    x = polynomial_ring.gen()
    q = x**2 + u*x + v
    H = h0 + h1*x + h2*x**2 + h3*x**3
    identity = H**2 - f - h3**2*q**3
    if identity[6] != 0:
        raise ArithmeticError("leading coefficient did not cancel")
    return tuple(identity[index] for index in range(6))


def symbolic_fiber_jacobian():
    """The fiber derivative is independent of the coefficients of f."""
    ring = PolynomialRing(ZZ, VARIABLE_NAMES)
    u, v, h0, h1, h2, h3 = ring.gens()
    polynomial_ring = PolynomialRing(ring, "X")
    X = polynomial_ring.gen()
    q = X**2 + u*X + v
    H = h0 + h1*X + h2*X**2 + h3*X**3
    identity = H**2 - h3**2*q**3
    equations = tuple(ring(identity[index]) for index in range(6))
    jacobian = matrix(
        ring,
        [[equation.derivative(variable) for variable in ring.gens()]
         for equation in equations],
    )
    return ring, equations, jacobian


SYMBOLIC_RING, SYMBOLIC_EQUATIONS, FIBER_JACOBIAN = symbolic_fiber_jacobian()


def evaluate_fiber_jacobian(field, values):
    substitutions = {
        variable: field(value)
        for variable, value in zip(SYMBOLIC_RING.gens(), values)
    }
    return FIBER_JACOBIAN.apply_map(
        lambda entry: field(entry.subs(substitutions))
    )


def enumerate_special_fiber():
    """Exhaust all normalized contact identities over F_7 at t=4."""
    field = GF(PRIME)
    _r, _a, _marked_h, f = contact9_polynomials(field, field(4))
    solutions = []
    for values in product(field, repeat=6):
        if all(value == 0 for value in fiber_equations(f, values)):
            solutions.append(tuple(Integer(value) for value in values))
    return f, solutions


def project_to_normalization(q, H):
    """Return the elliptic-normalization sum of the divisor defined by q,H."""
    field = GF(PRIME)
    extension = GF(PRIME**2, "z")
    polynomial_ring = PolynomialRing(extension, "x")
    x = polynomial_ring.gen()
    cubic = x**3 + 4*x**2 + 6*x + 5
    elliptic = EllipticCurve(extension, [0, 4, 0, 6, 5])
    q_extension = polynomial_ring(q)
    H_extension = polynomial_ring(H)
    points = [
        elliptic(alpha, H_extension(alpha)/(alpha + 1))
        for alpha, multiplicity in q_extension.roots(extension)
        for _unused in range(multiplicity)
    ]
    return sum(points, elliptic(0))


def generalized_ratio(q, oriented_H):
    """Compare D to 3*D9 in Pic^0 of the split nodal special curve.

    The orientation is chosen so the elliptic projection is Qminus=3*P.
    If psi has divisor D-(Qminus-O) and f3 is the Miller function with
    divisor 3(P-O)-(Qminus-O), then psi/f3 is principal for D-3D9.
    Equality in the generalized Jacobian is equivalent to equal values at
    the two points above the node.  The returned ratio is phi(Q+)/phi(Q-).
    """
    field = GF(PRIME)
    polynomial_ring = PolynomialRing(field, "x")
    x = polynomial_ring.gen()
    cubic = x**3 + 4*x**2 + 6*x + 5
    y_remainder = (oriented_H*(x + 1).inverse_mod(q)) % q
    intercept, slope = y_remainder[0], y_remainder[1]

    # P=(1,3), 2P=(3,3), 3P=Qminus=(6,4), and Qplus=(6,3).
    qminus = (field(6), field(4))
    qplus = (field(6), field(3))
    tangent_slope = field(4)
    second_line_slope = field(0)

    def candidate_line(point):
        xx, yy = point
        return yy - slope*xx - intercept

    def first_miller_line(point):
        xx, yy = point
        return yy - tangent_slope*xx - field(6)

    def second_miller_line(point):
        xx, yy = point
        return yy - field(3)

    phi_minus = (
        candidate_line(qminus)*(qminus[0] - field(3))
        / (first_miller_line(qminus)*second_miller_line(qminus))
    )

    # At Qplus both the candidate line and y-3 vanish.  Take their ratio
    # using x as a local parameter on the elliptic normalization.
    y_derivative = cubic.derivative()(qplus[0])/(2*qplus[1])
    cancelled_line_ratio = (
        (y_derivative - slope)/(y_derivative - second_line_slope)
    )
    phi_plus = (
        cancelled_line_ratio*(qplus[0] - field(3))
        / first_miller_line(qplus)
    )
    return phi_plus/phi_minus, slope, intercept


def classify_special_solutions(f, solutions):
    field = GF(PRIME)
    polynomial_ring = f.parent()
    x = polynomial_ring.gen()
    signed_pairs = {}
    for values in solutions:
        u, v, h0, h1, h2, h3 = map(field, values)
        q = x**2 + u*x + v
        H = h0 + h1*x + h2*x**2 + h3*x**3
        key = (Integer(u), Integer(v))
        signed_pairs.setdefault(key, []).append(H)

    records = []
    for key in sorted(signed_pairs):
        u, v = map(field, key)
        q = x**2 + u*x + v
        representatives = signed_pairs[key]
        H = representatives[0]
        jacobian = evaluate_fiber_jacobian(
            field, (u, v, H[0], H[1], H[2], H[3])
        )
        node_value = q(field(6))
        projection = None
        ratio = None
        oriented_H = None
        slope_intercept = None
        independent = None
        if node_value != 0:
            projection = project_to_normalization(q, H)
            projection_affine = (
                Integer(field(projection[0])), Integer(field(projection[1]))
            )
            if projection_affine == (6, 4):
                oriented_H = H
            elif projection_affine == (6, 3):
                oriented_H = -H
            else:
                raise ArithmeticError("unexpected elliptic projection")
            ratio, slope, intercept = generalized_ratio(q, oriented_H)
            slope_intercept = (Integer(slope), Integer(intercept))
            independent = ratio != 1
        records.append({
            "q_key": key,
            "q": q,
            "H": H,
            "signed_count": len(representatives),
            "node_value": Integer(node_value),
            "fiber_rank": jacobian.rank(),
            "fiber_det": Integer(jacobian.det()),
            "projection": projection,
            "ratio": None if ratio is None else Integer(ratio),
            "oriented_H": oriented_H,
            "line": slope_intercept,
            "independent": independent,
        })
    return records


def rational_t_values(height_bound):
    values = set()
    for denominator in range(1, height_bound + 1):
        if denominator % PRIME == 0:
            continue
        for numerator in range(-height_bound, height_bound + 1):
            if gcd(numerator, denominator) != 1:
                continue
            if (numerator - 4*denominator) % PRIME:
                continue
            value = QQ(numerator)/denominator
            if value not in (QQ(-1), QQ(1)):
                values.add(value)
    return sorted(values, key=lambda value: (
        max(abs(value.numerator()), value.denominator()),
        abs(value.numerator()) + value.denominator(),
        value,
    ))


def hensel_lift(base_t, residue, precision):
    """Uniquely lift one smooth fiber root modulo 7^precision."""
    field = GF(PRIME)
    inverse = evaluate_fiber_jacobian(field, residue).inverse()
    fiber = [Integer(value) for value in residue]
    modulus = Integer(PRIME)
    _r, _a, _marked_h, f = contact9_polynomials(QQ, QQ(base_t))

    for _exponent in range(1, precision):
        residuals = fiber_equations(f, tuple(QQ(value) for value in fiber))
        quotients = []
        for residual in residuals:
            residual = QQ(residual)
            numerator = Integer(residual.numerator())
            denominator = Integer(residual.denominator())
            if numerator % modulus:
                raise ArithmeticError("Hensel input is not a root at current precision")
            quotients.append(
                -field(numerator // modulus)/field(denominator)
            )
        correction = inverse*vector(field, quotients)
        next_modulus = modulus*PRIME
        fiber = [
            Integer((value + modulus*Integer(delta)) % next_modulus)
            for value, delta in zip(fiber, correction)
        ]
        modulus = next_modulus

    final_residuals = fiber_equations(f, tuple(QQ(value) for value in fiber))
    for residual in final_residuals:
        residual = QQ(residual)
        if Integer(residual.numerator()) % modulus:
            raise ArithmeticError("final Hensel residual is nonzero")
    return modulus, tuple(fiber)


def reconstruct_coordinate(residue, modulus):
    try:
        return QQ(Integer(residue).rational_reconstruction(modulus))
    except (ArithmeticError, ValueError):
        return None


def reconstruct_and_check(base_t, lift):
    modulus, residues = lift
    fiber = tuple(reconstruct_coordinate(value, modulus) for value in residues)
    if any(value is None for value in fiber):
        return fiber, None
    _r, _a, _marked_h, f = contact9_polynomials(QQ, QQ(base_t))
    residuals = fiber_equations(f, fiber)
    return fiber, residuals


def verify_exact_hit(base_t, fiber):
    r, a, marked_h, f = contact9_polynomials(QQ, QQ(base_t))
    polynomial_ring = f.parent()
    x = polynomial_ring.gen()
    u, v, h0, h1, h2, h3 = fiber
    q = x**2 + u*x + v
    H = h0 + h1*x + h2*x**2 + h3*x**3
    identity = H**2 - f - h3**2*q**3
    curve = HyperellipticCurve(f)
    jacobian = curve.jacobian()(QQ)
    zero = jacobian(0)
    candidate = jacobian([q, H % q])
    marked = jacobian([x - 1, polynomial_ring(marked_h(1))])
    root_two = jacobian([x - r, polynomial_ring(0)])
    relation = next(
        (coefficient for coefficient in range(9) if candidate == coefficient*marked),
        None,
    )
    return {
        "t": base_t,
        "a": a,
        "f": f,
        "q": q,
        "H": H,
        "identity": identity,
        "discriminant": f.discriminant(),
        "q_discriminant": q.discriminant(),
        "gcd_q_f": q.gcd(f),
        "candidate_order3": candidate != zero and 3*candidate == zero,
        "marked_order9": 9*marked == zero and 3*marked != zero,
        "root_order2": 2*root_two == zero and root_two != zero,
        "relation_to_marked": relation,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height-bound", type=int, default=80)
    parser.add_argument("--precision", type=int, default=14)
    parser.add_argument("--max-t", type=int, default=0,
                        help="0 means every t in the declared height box")
    args = parser.parse_args()

    print("# contact-9 independent order-3 cover at t=4 mod 7")
    print("identity H^2-f=h3^2*q^3, q=x^2+u*x+v, deg(H)=3")
    print("six coefficient equations in six fiber variables over the t-line")
    f_special, solutions = enumerate_special_fiber()
    print(f"special_f={f_special}")
    print(f"special_factorization={f_special.factor()}")
    print(f"signed_solution_count={len(solutions)}")
    print(f"signed_solutions={solutions}")

    records = classify_special_solutions(f_special, solutions)
    print("special_q_records")
    for record in records:
        print(
            " q=%s H=%s signed=%s q_at_node=%s fiber_rank=%s fiber_det=%s "
            "projection=%s oriented_H=%s line=%s generalized_ratio=%s independent=%s"
            % (
                record["q"], record["H"], record["signed_count"],
                record["node_value"], record["fiber_rank"], record["fiber_det"],
                record["projection"], record["oriented_H"], record["line"],
                record["ratio"], record["independent"],
            )
        )

    independent_records = [record for record in records if record["independent"]]
    if len(independent_records) != 2:
        raise ArithmeticError("expected exactly two independent signed q-pairs")
    branch_residues = [
        tuple(Integer(value) for value in (
            record["q"][1], record["q"][0], record["H"][0],
            record["H"][1], record["H"][2], record["H"][3],
        ))
        for record in independent_records
    ]
    print(f"independent_branch_residues={branch_residues}")

    base_values = rational_t_values(args.height_bound)
    if args.max_t:
        base_values = base_values[:args.max_t]
    print(f"height_bound={args.height_bound}")
    print(f"precision={args.precision}")
    print(f"modulus={PRIME**args.precision}")
    reconstruction_bound = isqrt((PRIME**args.precision)//2)
    print(f"componentwise_reconstruction_height_bound={reconstruction_bound}")
    print(f"t_count={len(base_values)}")
    print(f"first_t_values={base_values[:20]}")

    exact_hits = []
    branch_counts = Counter()
    sample_lifts = []
    for t_index, base_t in enumerate(base_values):
        for branch_index, residue in enumerate(branch_residues):
            lift = hensel_lift(base_t, residue, args.precision)
            fiber, residuals = reconstruct_and_check(base_t, lift)
            if base_t == 4:
                sample_lifts.append((branch_index, lift, fiber, residuals))
            if residuals is None:
                branch_counts[(branch_index, "incomplete_reconstruction")] += 1
                continue
            branch_counts[(branch_index, "complete_reconstruction")] += 1
            if any(residuals):
                continue
            verification = verify_exact_hit(base_t, fiber)
            exact_hits.append((t_index, branch_index, fiber, verification))
            print(f"EXACT_HIT t_index={t_index} branch_index={branch_index}")
            print(f"fiber={fiber}")
            for key, value in verification.items():
                print(f" {key}={value}")

    print("t_equals_4_lift_samples")
    for branch_index, lift, fiber, residuals in sample_lifts:
        print(
            f" branch={branch_index} residues_mod_{lift[0]}={lift[1]} "
            f"reconstruction={fiber} exact_residuals={residuals}"
        )
    print(f"branch_reconstruction_counts={dict(sorted(branch_counts.items()))}")
    print(f"tested_lifts={len(base_values)*len(branch_residues)}")
    print(f"exact_hit_count={len(exact_hits)}")
    if not exact_hits:
        print("BOUNDED_RESULT=NO_RATIONAL_POINT_RECONSTRUCTED_ON_TWO_INDEPENDENT_BRANCHES")


if __name__ == "__main__":
    main()
