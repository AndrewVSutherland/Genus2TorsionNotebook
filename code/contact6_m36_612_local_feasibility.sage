#!/usr/bin/env sage
"""Local feasibility of the contact-6 route to Z/6 x Z/12.

Intersect the cubic-contact equations F1=F2=F3 with the two equations
obtained from H1=...=H4 after eliminating the linear variables in the
marked-class halving identity.  At small primes the script

* enumerates the complete affine intersection;
* removes the explicit degeneracy locus pointwise;
* requires the [1,2,2]-compatible factorisation used by the known [6,6]
  example;
* recovers the full finite Jacobian invariants, rather than screening only
  by its order;
* verifies that the displayed contact and halving classes generate a
  subgroup Z/6 x Z/12; and
* counts all lifts modulo p^2 by linear Hensel equations.

Typical run:
    sage code/contact6_m36_612_local_feasibility.sage --primes 5,7
"""

import argparse
from collections import Counter, defaultdict
from itertools import combinations, product

from sage.all import (
    GF,
    HyperellipticCurve,
    Matrix,
    PolynomialRing,
    QQ,
    ZZ,
    matrix,
    prime_divisors,
    prod,
    vector,
)


TARGET = [ZZ(6), ZZ(12)]


def strip_p_part(n, p):
    n = ZZ(n)
    while n % p == 0:
        n //= p
    return n


def prime_to_p_target(p):
    return [m for n in TARGET if (m := strip_p_part(n, p)) > 1]


def has_subgroup_embedding(ambient, required):
    """Test finite abelian subgroup containment prime by prime."""
    if not required:
        return True
    for ell in prime_divisors(prod(required)):
        max_exp = max(ZZ(n).valuation(ell) for n in required)
        for k in range(1, max_exp + 1):
            need = sum(ZZ(n).valuation(ell) >= k for n in required)
            have = sum(ZZ(n).valuation(ell) >= k for n in ambient)
            if have < need:
                return False
    return True


def square_roots(value):
    if not value.is_square():
        return []
    if value == 0:
        return [value]
    root = value.sqrt()
    return [root, -root]


def reduced_jacobian_points(f):
    """Enumerate the unique reduced Mumford pairs for y^2=f, deg(f)=5."""
    field = f.base_ring()
    ring = f.parent()
    x = ring.gen()
    jacobian_points = HyperellipticCurve(f).jacobian()(field)
    points_by_rep = {repr(jacobian_points(0)): jacobian_points(0)}

    def add_point(divisor):
        points_by_rep[repr(divisor)] = divisor

    for a in field:
        u = x + a
        for d in square_roots(f(-a)):
            add_point(jacobian_points([u, ring(d)]))

    aux = PolynomialRing(field, "w")
    w = aux.gen()
    for a in field:
        for b in field:
            u = x**2 + a*x + b
            remainder = f % u
            coeff_x = remainder[1]
            constant = remainder[0]
            candidates = set()
            if coeff_x == 0:
                for d in square_roots(constant):
                    candidates.add((field(0), d))

            equation = (
                (a**2 - 4*b)*w**2
                + (2*a*coeff_x - 4*constant)*w
                + coeff_x**2
            )
            if equation:
                for value, _multiplicity in equation.roots():
                    if value == 0:
                        continue
                    for c in square_roots(value):
                        d = (coeff_x + a*c**2)/(2*c)
                        candidates.add((c, d))

            for c, d in candidates:
                vpoly = c*x + d
                if (f - vpoly**2) % u:
                    raise ArithmeticError("Mumford solution failed verification")
                try:
                    add_point(jacobian_points([u, vpoly]))
                except (ArithmeticError, ValueError, ZeroDivisionError):
                    # A non-coprime pair represents a divisor of lower degree.
                    pass
    return jacobian_points, list(points_by_rep.values())


def jacobian_order_from_curve(f):
    curve = HyperellipticCurve(f)
    n1, n2 = curve.count_points(2)
    p = f.base_ring().cardinality()
    s1 = p + 1 - n1
    sum_alpha2 = p**2 + 1 - n2
    e2 = (s1**2 - sum_alpha2)//2
    return ZZ(1 - s1 + e2 - p*s1 + p**2)


def abelian_invariants_and_points(f):
    jacobian_points, points = reduced_jacobian_points(f)
    expected = jacobian_order_from_curve(f)
    if len(points) != expected:
        raise ArithmeticError(
            f"Mumford enumeration has {len(points)} points, expected {expected}"
        )

    zero = jacobian_points(0)
    sylow_exponents = {}
    for ell, total_exp in expected.factor():
        kernel_logs = []
        for k in range(1, total_exp + 1):
            kernel_size = sum(ell**k*divisor == zero for divisor in points)
            log_size = ZZ(kernel_size).valuation(ell)
            if ell**log_size != kernel_size:
                raise ArithmeticError("multiplication kernel is not an ell-power")
            kernel_logs.append(log_size)
            if log_size == total_exp:
                break

        ranks = []
        previous = 0
        for log_size in kernel_logs:
            ranks.append(log_size - previous)
            previous = log_size
        ranks.append(0)

        exponents = []
        for k in range(1, len(ranks)):
            exponents.extend([k]*(ranks[k - 1] - ranks[k]))
        sylow_exponents[ZZ(ell)] = sorted(exponents)

    rank = max((len(exponents) for exponents in sylow_exponents.values()), default=0)
    invariants = [ZZ(1)]*rank
    for ell, exponents in sylow_exponents.items():
        padded = [0]*(rank - len(exponents)) + exponents
        invariants = [n*ell**exponent for n, exponent in zip(invariants, padded)]
    return [n for n in invariants if n > 1], jacobian_points, points


def point_order(point, group_order):
    if not point:
        return ZZ(1)
    order = ZZ(group_order)
    for ell, exponent in order.factor():
        for _ in range(exponent):
            if (order//ell)*point == 0:
                order //= ell
            else:
                break
    return order


def contact_and_halving_data(field, values):
    """Return the five equations and all curve/divisor data."""
    a, b, contact_scale, contact_u, contact_v, half_q, half_r = values
    ring = PolynomialRing(field, "x")
    x = ring.gen()
    t = x - 1

    h6 = 1 + a*x + b*x**2 + x**3
    f = h6**2 - (x - 1)**6
    c1, c2, c3, c4, c5 = [f[i] for i in range(1, 6)]
    scale_square = contact_scale**2
    contact_b = c5*scale_square + 3*contact_u
    delta = (
        4*c4*scale_square + 12*(contact_u**2 + contact_v**2)
        - contact_b**2
    )
    f3 = (
        contact_b*delta + 16*contact_v**3 - 8*c3*scale_square
        - 8*contact_u**3 - 48*contact_u*contact_v**2
    )
    f2 = (
        delta**2 + 64*contact_b*contact_v**3 - 64*c2*scale_square
        - 192*(contact_u**2*contact_v**2 + contact_v**4)
    )
    f1 = delta*contact_v**3 - 4*c1*scale_square - 12*contact_u*contact_v**4

    center_a = b + 3
    center_c = a + b + 2
    k2 = (
        center_a**2*center_c**2*half_q
        - 2*center_a**2*center_c*half_q**2
        + center_a**2*half_q**4
        - 2*center_a*center_c**3
        - 2*center_a*center_c**2*half_q**2
        + 2*center_a*center_c**2*half_q
        + 4*center_a*center_c*half_q**2
        + 2*center_c**3*half_q
        - 2*center_c**3*half_r
        - center_c**2*half_q*half_r**2
        - 4*center_c**2*half_q
    )
    # Eight times K3 from contact6_m36_halveD_symbolic.m.
    k3x8 = (
        center_a**4*center_c - 12*center_a**3*center_c
        - 12*center_a**2*center_c**2 + 32*center_a**2*center_c*half_q
        - 2*center_a**2*center_c*half_r**2
        - 16*center_a**2*center_c*half_r + 28*center_a**2*center_c
        + 16*center_a**2*half_q**2*half_r
        - 16*center_a*center_c**2*half_r - 8*center_a*center_c**2
        - 4*center_a*center_c*half_r**2 + 32*center_a*center_c*half_r
        - 16*center_a*center_c + 4*center_c**3
        - 4*center_c**2*half_r**2 - 16*center_c**2
        + center_c*half_r**4 + 8*center_c*half_r**2 + 16*center_c
    )

    data = {
        "ring": ring,
        "x": x,
        "t": t,
        "h6": h6,
        "f": f,
        "contact_b": contact_b,
        "delta": delta,
        "contact_q": x**2 + contact_u*x + contact_v**2,
        "center_a": center_a,
        "center_c": center_c,
    }
    return (f1, f2, f3, k2, k3x8), data


def boundary_labels(data, values):
    a, b, contact_scale, contact_u, contact_v, half_q, _half_r = values
    labels = []
    if data["center_a"] == 0:
        labels.append("A0")
    if data["center_c"] == 0:
        labels.append("C0")
    if contact_scale == 0:
        labels.append("L0")
    if contact_v == 0:
        labels.append("v0")
    if contact_v**3 == 1:
        labels.append("v3eq1")
    if contact_u**2 == 4*contact_v**2:
        labels.append("contact_q_double")
    if half_q == 0:
        labels.append("half_q0")
    f = data["f"]
    if f.degree() != 5:
        labels.append("degree_drop")
    elif f.discriminant() == 0:
        labels.append("singular_curve")
    if data["contact_q"].gcd(f).degree() > 0:
        labels.append("contact_q_meets_f")
    return labels


def recover_half(data, values):
    field = data["f"].base_ring()
    _a, _b, _contact_scale, _contact_u, _contact_v, half_q, half_r = values
    center_a = data["center_a"]
    center_c = data["center_c"]
    if center_a == 0 or center_c == 0:
        return None
    center_b = center_a + center_c - 2
    half_p = (center_a**2 + 2*center_b - half_r**2)/(4*center_a)
    half_s = center_a*half_q**2/center_c - center_b
    t = data["t"]
    half_u = t**2 + half_p*t + half_q
    half_ell = half_r*t**2 + half_s*t - center_c
    identity = (
        half_ell**2 - data["f"]
        + 2*center_a*data["t"]*half_u**2
    )
    if identity:
        raise ArithmeticError("eliminated halving equations failed verification")
    return half_p, half_s, half_u, half_ell


def core_factor_compatible(f):
    if f.degree() != 5 or f.discriminant() == 0:
        return False
    degrees = sorted(factor.degree() for factor, _exponent in f.factor())
    # A Q-factorisation x*(quadratic)*(quadratic) can only reduce with
    # irreducible degrees 1 or 2 at a good prime.
    return len(degrees) >= 3 and max(degrees) <= 2


def geometrically_nonsplit_certificate(f):
    curve = HyperellipticCurve(f)
    frobenius = curve.frobenius_polynomial()
    bivariate = PolynomialRing(QQ, 2, "zw")
    z, w = bivariate.gens()
    frobenius_bivariate = sum(QQ(frobenius[i])*z**i for i in range(5))
    resultant = frobenius_bivariate.resultant(w - z**12, z)
    power_polynomial = PolynomialRing(QQ, "w")(resultant)
    return power_polynomial.is_irreducible(), frobenius, power_polynomial


def two_torsion_points(jacobian_points, f):
    ring = f.parent()
    factors = [factor.monic() for factor, _exponent in f.factor()]
    points = []
    seen = set()
    for size in range(1, len(factors) + 1):
        for subset in combinations(factors, size):
            u = prod(subset, ring(1)).monic()
            if u.degree() > 2:
                continue
            try:
                divisor = jacobian_points([u, ring(0)])
            except (ArithmeticError, ValueError, ZeroDivisionError):
                continue
            key = repr(divisor)
            if key not in seen and divisor != 0:
                seen.add(key)
                points.append(divisor)
    return points


def explicit_target_check(data, values, invariants, jacobian_points):
    field = data["f"].base_ring()
    ring = data["ring"]
    x = data["x"]
    a, b, contact_scale, contact_u, contact_v, _half_q, _half_r = values
    recovered = recover_half(data, values)
    if recovered is None:
        return False, None
    _half_p, _half_s, half_u, half_ell = recovered
    if half_u.gcd(data["f"]).degree() > 0:
        return False, None

    contact_q = data["contact_q"]
    if contact_scale == 0:
        return False, None
    h3 = (
        x**3 + data["contact_b"]/2*x**2 + data["delta"]/8*x
        + contact_v**3
    )/contact_scale
    try:
        marked_d = jacobian_points([x - 1, ring(data["center_c"])])
        contact_e = jacobian_points([contact_q, h3 % contact_q])
        half_h = jacobian_points([half_u, half_ell % half_u])
    except (ArithmeticError, ValueError, ZeroDivisionError):
        return False, None

    group_order = prod(invariants, ZZ(1))
    ord_d = point_order(marked_d, group_order)
    ord_e = point_order(contact_e, group_order)
    ord_h = point_order(half_h, group_order)
    double_sign = 1 if 2*half_h == marked_d else (-1 if 2*half_h == -marked_d else 0)
    if (ord_d, ord_e, ord_h, abs(double_sign)) != (6, 3, 12, 1):
        return False, (ord_d, ord_e, ord_h, double_sign, 0)

    for torsion_2 in two_torsion_points(jacobian_points, data["f"]):
        second = contact_e + torsion_2
        if point_order(second, group_order) != 6:
            continue
        generated = {
            repr(i*half_h + j*second)
            for i in range(12)
            for j in range(6)
        }
        if len(generated) == 72:
            return True, (ord_d, ord_e, ord_h, double_sign, len(generated))
    return False, (ord_d, ord_e, ord_h, double_sign, 0)


def integer_equations():
    ring = PolynomialRing(ZZ, ("a", "b", "L", "U", "v", "q", "r"))
    a, b, contact_scale, contact_u, contact_v, half_q, half_r = ring.gens()
    xring = PolynomialRing(ring, "x")
    x = xring.gen()
    h6 = 1 + a*x + b*x**2 + x**3
    f = h6**2 - (x - 1)**6
    c1, c2, c3, c4, c5 = [f[i] for i in range(1, 6)]
    scale_square = contact_scale**2
    contact_b = c5*scale_square + 3*contact_u
    delta = 4*c4*scale_square + 12*(contact_u**2 + contact_v**2) - contact_b**2
    f3 = contact_b*delta + 16*contact_v**3 - 8*c3*scale_square - 8*contact_u**3 - 48*contact_u*contact_v**2
    f2 = delta**2 + 64*contact_b*contact_v**3 - 64*c2*scale_square - 192*(contact_u**2*contact_v**2 + contact_v**4)
    f1 = delta*contact_v**3 - 4*c1*scale_square - 12*contact_u*contact_v**4
    center_a = b + 3
    center_c = a + b + 2
    k2 = center_a**2*center_c**2*half_q - 2*center_a**2*center_c*half_q**2 + center_a**2*half_q**4 - 2*center_a*center_c**3 - 2*center_a*center_c**2*half_q**2 + 2*center_a*center_c**2*half_q + 4*center_a*center_c*half_q**2 + 2*center_c**3*half_q - 2*center_c**3*half_r - center_c**2*half_q*half_r**2 - 4*center_c**2*half_q
    k3x8 = center_a**4*center_c - 12*center_a**3*center_c - 12*center_a**2*center_c**2 + 32*center_a**2*center_c*half_q - 2*center_a**2*center_c*half_r**2 - 16*center_a**2*center_c*half_r + 28*center_a**2*center_c + 16*center_a**2*half_q**2*half_r - 16*center_a*center_c**2*half_r - 8*center_a*center_c**2 - 4*center_a*center_c*half_r**2 + 32*center_a*center_c*half_r - 16*center_a*center_c + 4*center_c**3 - 4*center_c**2*half_r**2 - 16*center_c**2 + center_c*half_r**4 + 8*center_c*half_r**2 + 16*center_c
    equations = [ring(poly) for poly in (f1, f2, f3, k2, k3x8)]
    return ring, equations


INTEGER_RING, INTEGER_EQUATIONS = integer_equations()
INTEGER_JACOBIAN = Matrix(INTEGER_RING, [
    [equation.derivative(variable) for variable in INTEGER_RING.gens()]
    for equation in INTEGER_EQUATIONS
])


def hensel_information(p, residue):
    """Count and exhibit all linear lifts of a residue point modulo p^2."""
    field = GF(p)
    integer_point = [ZZ(value) for value in residue]
    constants = []
    for equation in INTEGER_EQUATIONS:
        value = ZZ(equation(*integer_point))
        if value % p:
            raise ArithmeticError("input is not a solution modulo p")
        constants.append(field(-(value//p)))
    jacobian = matrix(
        field,
        [[entry(*integer_point) for entry in row] for row in INTEGER_JACOBIAN.rows()],
    )
    rhs = vector(field, constants)
    augmented = jacobian.augment(rhs.column())
    rank = jacobian.rank()
    if augmented.rank() != rank:
        return rank, 0, None
    tangent = jacobian.right_kernel()
    lift_count = ZZ(p)**tangent.dimension()
    correction = jacobian.solve_right(rhs)
    modulus = p**2
    lift = tuple((integer_point[i] + p*ZZ(correction[i])) % modulus for i in range(7))
    if any(ZZ(equation(*lift)) % modulus for equation in INTEGER_EQUATIONS):
        raise ArithmeticError("reported Hensel lift failed verification")
    return rank, lift_count, lift


def finite_study(p, show_limit):
    field = GF(p)
    required = prime_to_p_target(p)
    print(f"\n=== prime {p}; required prime-to-p subgroup {required} ===")

    equation_solutions = []
    boundary_counts = Counter()
    open_points = []
    smooth_core_points = []
    cache = {}
    ambient_counts = Counter()
    order_passes = embedding_passes = explicit_passes = 0
    nonsplit_passes = 0
    false_order_positives = 0
    target_records = []
    target_frobenius_counts = Counter()
    target_power12_factor_counts = Counter()
    target_lift_counts = Counter()

    # Complete affine enumeration.  Filtering F1,F2,F3 before q,r keeps the
    # seven-variable loop small even at p=7.
    for first_five in product(field, repeat=5):
        seed = tuple(first_five) + (field(0), field(0))
        contact_equations, _data = contact_and_halving_data(field, seed)
        if any(contact_equations[i] for i in range(3)):
            continue
        for half_q, half_r in product(field, repeat=2):
            values = tuple(first_five) + (half_q, half_r)
            equations, data = contact_and_halving_data(field, values)
            if equations[3] or equations[4]:
                continue
            equation_solutions.append(values)
            labels = boundary_labels(data, values)
            signature = "+".join(sorted(labels)) if labels else "equation_open"
            boundary_counts[signature] += 1

            f = data["f"]
            if f.degree() == 5 and f.discriminant() != 0 and core_factor_compatible(f):
                nonsplit, frobenius, _power = geometrically_nonsplit_certificate(f)
                if nonsplit:
                    smooth_core_points.append((values, data, labels, frobenius))

            if labels or not core_factor_compatible(f):
                continue
            recovered = recover_half(data, values)
            if recovered is None or recovered[2].gcd(f).degree() > 0:
                continue
            open_points.append((values, data))

            key = tuple(ZZ(f[i]) for i in range(6))
            if key not in cache:
                cache[key] = abelian_invariants_and_points(f)
            invariants, jacobian_points, _points = cache[key]
            ambient_counts[tuple(invariants)] += 1
            required_order = prod(required, ZZ(1))
            group_order = prod(invariants, ZZ(1))
            order_ok = group_order % required_order == 0
            embedding_ok = has_subgroup_embedding(invariants, required)
            order_passes += order_ok
            embedding_passes += embedding_ok
            false_order_positives += order_ok and not embedding_ok
            if not embedding_ok:
                continue
            explicit_ok, class_data = explicit_target_check(
                data, values, invariants, jacobian_points
            )
            if not explicit_ok:
                continue
            explicit_passes += 1
            nonsplit, frobenius, power_polynomial = geometrically_nonsplit_certificate(f)
            nonsplit_passes += nonsplit
            target_frobenius_counts[str(frobenius.factor())] += 1
            target_power12_factor_counts[str(power_polynomial.factor())] += 1
            target_hensel = hensel_information(p, values)
            target_lift_counts[(target_hensel[0], target_hensel[1])] += 1
            if nonsplit:
                record = {
                    "values": values,
                    "invariants": invariants,
                    "classes": class_data,
                    "frobenius": frobenius,
                    "power": power_polynomial,
                }
                record["hensel"] = hensel_information(p, values)
                target_records.append(record)
    component_counts = Counter()
    for signature, count in boundary_counts.items():
        if signature == "equation_open":
            component_counts[signature] += count
            continue
        for label in signature.split("+"):
            component_counts[label] += count
    top_signatures = sorted(
        boundary_counts.items(), key=lambda item: (-item[1], item[0])
    )[:12]
    print(f"affine_intersection_points {len(equation_solutions)}")
    print(f"boundary_component_incidences {dict(sorted(component_counts.items()))}")
    print(f"largest_boundary_signatures {top_signatures}")
    print(f"smooth_core_nonsplit_points {len(smooth_core_points)}")
    print(f"fully_open_core_points {len(open_points)}")
    print(f"ambient_invariant_counts {dict(sorted(ambient_counts.items()))}")
    print(
        f"order_passes {order_passes} embedding_passes {embedding_passes} "
        f"false_order_positives {false_order_positives} "
        f"explicit_6x12_passes {explicit_passes} "
        f"geometrically_nonsplit_passes {nonsplit_passes}"
    )

    print(f"target_frobenius_factor_counts {dict(sorted(target_frobenius_counts.items()))}")
    print(f"target_power12_factor_counts {dict(sorted(target_power12_factor_counts.items()))}")
    print(f"target_hensel_rank_lift_counts {dict(sorted(target_lift_counts.items()))}")

    rank_counts = Counter()
    lift_counts = Counter()
    smooth_core_lift_counts = Counter()
    smooth_core_keys = {tuple(point[0]) for point in smooth_core_points}
    samples = []
    for values in equation_solutions:
        rank, count, lift = hensel_information(p, values)
        rank_counts[rank] += 1
        lift_counts[count] += 1
        if tuple(values) in smooth_core_keys:
            smooth_core_lift_counts[(rank, count)] += 1
        if len(samples) < show_limit and count:
            samples.append((tuple(ZZ(v) for v in values), rank, count, lift))
    print(f"intersection_jacobian_rank_counts {dict(sorted(rank_counts.items()))}")
    print(f"per_point_lift_count_distribution {dict(sorted(lift_counts.items()))}")
    print(f"smooth_core_nonsplit_lift_distribution {dict(sorted(smooth_core_lift_counts.items()))}")
    print(f"intersection_lift_samples {samples}")

    for index, record in enumerate(target_records[:show_limit], start=1):
        print(
            f"TARGET_NONSPLIT {index} values {tuple(ZZ(v) for v in record['values'])} "
            f"invariants {record['invariants']} classes {record['classes']} "
            f"frobenius {record['frobenius']} hensel {record['hensel']}"
        )
    return {
        "intersection": len(equation_solutions),
        "smooth_core_nonsplit": len(smooth_core_points),
        "open": len(open_points),
        "embedding": embedding_passes,
        "explicit": explicit_passes,
        "target_nonsplit": len(target_records),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="5,7")
    parser.add_argument("--show-limit", type=int, default=8)
    args = parser.parse_args()
    primes = [ZZ(part.strip()) for part in args.primes.split(",") if part.strip()]
    print("Contact-6 [6,6] core plus marked-class halving local feasibility")
    print(f"target {TARGET} primes {primes}")
    print("variables (a,b,L,U,v,q,r); equations (F1,F2,F3,K2,8*K3)")
    print("core condition: good degree-5 curve and [1,2,2]-compatible factor degrees")
    summaries = {p: finite_study(p, args.show_limit) for p in primes}
    print(f"\nSUMMARY {summaries}")


if __name__ == "__main__":
    main()
