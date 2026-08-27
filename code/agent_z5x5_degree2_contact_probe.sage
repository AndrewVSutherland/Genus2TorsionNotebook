#!/usr/bin/env sage
"""Bounded degree-2 contact probe for two independent 5-torsion classes.

This is intentionally a finite-field geometry test, not a rational search.
For

    f = (1 + a*x + b*x^2)^2 - k*x^5

the marked point (0,1) gives a rational class D0 of order 5.  A second
degree-2 class D=(q,w), q=x^2+u*x+v, has order 5 when a function
H-R*y has divisor 5D.  After normalising H to be monic this is the
polynomial Pell/contact identity

    H^2 - f*R^2 = q^5,

with deg(H)=5 and deg(R)<=2.

The script finds one open solution at each requested prime, verifies the
Jacobian classes, computes the equation-scheme tangent rank, and applies
the 12th-power Frobenius irreducibility test.  All loops are bounded by the
finite fields in PRIMES.
"""

from itertools import product

from sage.all import GF, HyperellipticCurve, PolynomialRing, QQ, matrix


PRIMES = (11, 19)


def power_frobenius_certificate(curve):
    """Return Frobenius data and the standard 12th-power transform test."""
    frobenius = curve.frobenius_polynomial()
    bivariate = PolynomialRing(QQ, 2, "zw")
    z, w = bivariate.gens()
    frobenius_bivariate = sum(QQ(frobenius[i])*z**i for i in range(5))
    resultant = frobenius_bivariate.resultant(w - z**12, z)
    power_polynomial = PolynomialRing(QQ, "w")(resultant)
    return frobenius, power_polynomial, power_polynomial.is_irreducible()


def monic_square_root(poly, x):
    """Return the monic degree-5 square root of poly, or None."""
    if poly.degree() != 10 or poly[10] != 1:
        return None
    field = poly.base_ring()
    root = x**5
    for degree in range(9, 4, -1):
        coefficient = (poly[degree] - (root**2)[degree]) / field(2)
        root += coefficient*x**(degree - 5)
    return root if root**2 == poly else None


def relation_coefficients(point, marked):
    """Return c in F_5 with point=c*marked, or None if independent."""
    for coefficient in range(5):
        if point == coefficient*marked:
            return coefficient
    return None


def recover_contact_identity(f, q, w, x):
    """Boundedly scan R and recover H^2-f*R^2=q^5."""
    field = f.base_ring()
    for r2, r1, r0 in product(field, repeat=3):
        Rpoly = r0 + r1*x + r2*x**2
        if Rpoly.gcd(q).degree() > 0:
            continue
        Hpoly = monic_square_root(q**5 + f*Rpoly**2, x)
        if Hpoly is None:
            continue
        minus_contact = (Hpoly - Rpoly*w) % q == 0
        plus_contact = (Hpoly + Rpoly*w) % q == 0
        if minus_contact or plus_contact:
            return Hpoly, Rpoly, "D" if minus_contact else "-D"
    return None


def second_class(curve, f, marked, x):
    """Find an independent order-5 class with irreducible squarefree q."""
    field = f.base_ring()
    jacobian_points = curve.jacobian()(field)
    zero = jacobian_points(0)
    for u, v in product(field, repeat=2):
        q = x**2 + u*x + v
        if q.is_irreducible() is False:
            continue
        if q.gcd(f).degree() > 0:
            continue
        for w1, w0 in product(field, repeat=2):
            w = w0 + w1*x
            if (w**2 - f) % q != 0:
                continue
            try:
                point = jacobian_points([q, w])
            except (ArithmeticError, TypeError, ValueError, ZeroDivisionError):
                continue
            if point == zero or 5*point != zero:
                continue
            if relation_coefficients(point, marked) is not None:
                continue
            recovered = recover_contact_identity(f, q, w, x)
            if recovered is not None:
                Hpoly, Rpoly, contact_sign = recovered
                return q, w, point, Hpoly, Rpoly, contact_sign
    return None


def equation_rank(prime, values):
    """Evaluate the ten coefficient equations and their Jacobian ranks."""
    field = GF(prime)
    names = (
        "a", "b", "k", "u", "v", "r0", "r1", "r2",
        "s0", "s1", "s2", "s3", "s4",
    )
    ring = PolynomialRing(field, names, order="degrevlex")
    (a, b, k, u, v, r0, r1, r2, s0, s1, s2, s3, s4) = ring.gens()
    poly_ring = PolynomialRing(ring, "X")
    X = poly_ring.gen()
    f = (1 + a*X + b*X**2)**2 - k*X**5
    q = X**2 + u*X + v
    Rpoly = r0 + r1*X + r2*X**2
    Hpoly = X**5 + s4*X**4 + s3*X**3 + s2*X**2 + s1*X + s0
    identity = Hpoly**2 - f*Rpoly**2 - q**5
    equations = [ring(identity[i]) for i in range(10)]
    point = {generator: field(value) for generator, value in zip(ring.gens(), values)}
    residuals = [equation.subs(point) for equation in equations]
    jacobian = matrix(
        ring,
        [[equation.derivative(variable) for variable in ring.gens()]
         for equation in equations],
    )
    evaluated = jacobian.apply_map(lambda entry: field(entry.subs(point)))
    fiber = evaluated.matrix_from_columns(range(3, 13))
    return residuals, evaluated.rank(), fiber.rank(), fiber.det(), equations


def find_open_solution(prime):
    field = GF(prime)
    polynomial_ring = PolynomialRing(field, "x")
    x = polynomial_ring.gen()

    tested = 0
    smooth = 0
    count_25 = 0
    absolute_simple = 0
    for a, b, k in product(field, repeat=3):
        if k == 0:
            continue
        tested += 1
        h0 = 1 + a*x + b*x**2
        f = h0**2 - k*x**5
        if f.degree() != 5 or f.discriminant() == 0:
            continue
        smooth += 1
        curve = HyperellipticCurve(f)
        frobenius, power_polynomial, is_absolute_simple = power_frobenius_certificate(curve)
        jacobian_order = int(frobenius(1))
        if jacobian_order % 25 != 0:
            continue
        count_25 += 1
        if not is_absolute_simple:
            continue
        absolute_simple += 1
        jacobian_points = curve.jacobian()(field)
        marked = jacobian_points([x, polynomial_ring(1)])
        if 5*marked != jacobian_points(0):
            raise ArithmeticError("marked contact class failed its order-5 identity")
        found = second_class(curve, f, marked, x)
        if found is None:
            continue
        q, w, point, Hpoly, Rpoly, contact_sign = found

        values = [
            a, b, k, q[1], q[0], Rpoly[0], Rpoly[1], Rpoly[2],
            Hpoly[0], Hpoly[1], Hpoly[2], Hpoly[3], Hpoly[4],
        ]
        residuals, full_rank, fiber_rank, fiber_det, equations = equation_rank(
            prime, values
        )
        if any(residuals):
            raise ArithmeticError("symbolic coefficient equations failed")
        return {
            "prime": prime,
            "tested": tested,
            "smooth": smooth,
            "count_25": count_25,
            "absolute_simple": absolute_simple,
            "a": a,
            "b": b,
            "k": k,
            "f": f,
            "q": q,
            "w": w,
            "H": Hpoly,
            "R": Rpoly,
            "contact_sign": contact_sign,
            "marked": marked,
            "second": point,
            "relation": relation_coefficients(point, marked),
            "frobenius": frobenius,
            "jacobian_order": jacobian_order,
            "power_polynomial": power_polynomial,
            "absolute_simple_certificate": is_absolute_simple,
            "full_rank": full_rank,
            "tangent_dimension": 13 - full_rank,
            "fiber_rank": fiber_rank,
            "fiber_det": fiber_det,
            "equations": equations,
        }
    return None


print("# degree-2 order-5 contact equations")
print("# f=(1+a*x+b*x^2)^2-k*x^5")
print("# q=x^2+u*x+v, R=r0+r1*x+r2*x^2")
print("# H=x^5+s4*x^4+s3*x^3+s2*x^2+s1*x+s0")
print("# E_i = coeff_x^i(H^2-f*R^2-q^5), i=0,...,9")
print("# There are 13 variables and 10 equations; expected dimension is 3.")

all_results = []
for prime in PRIMES:
    result = find_open_solution(prime)
    if result is None:
        print(f"p={prime}: NO_OPEN_SOLUTION_IN_BOUNDED_SCAN")
        continue
    all_results.append(result)
    print(f"\n## p={prime}")
    for key in (
        "tested", "smooth", "count_25", "absolute_simple", "a", "b", "k",
        "f", "q", "w", "H", "R", "contact_sign", "marked", "second",
        "relation", "frobenius", "jacobian_order", "power_polynomial",
        "absolute_simple_certificate", "full_rank", "tangent_dimension",
        "fiber_rank", "fiber_det",
    ):
        print(f"{key}={result[key]}")
    print(f"q_discriminant={result['q'].discriminant()}")
    print(f"q_irreducible={result['q'].is_irreducible()}")
    print(f"gcd_q_f={result['q'].gcd(result['f'])}")
    print(f"gcd_q_R={result['q'].gcd(result['R'])}")
    print(f"identity_check={result['H']**2-result['f']*result['R']**2-result['q']**5}")

if all_results:
    print("\n# expanded coefficient equations (over the final residue field)")
    for index, equation in enumerate(all_results[-1]["equations"]):
        print(f"E{index}={equation}")
