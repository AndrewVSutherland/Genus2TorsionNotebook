#!/usr/bin/env sage
"""Bounded exact-fiber and two-prime lifting test for contact-6 [2,6,6].

This is the executable Sage/Singular counterpart of
``contact6_m36_266_targeted_lift.m``.  It deliberately studies fixed fibers,
not a growing rational-height box.
"""

from argparse import ArgumentParser
from itertools import product
from math import gcd, isqrt

from sage.all import GF, QQ, ZZ, PolynomialRing, crt, lcm, matrix, vector


def primitive_integral(g):
    den = lcm([QQ(c).denominator() for c in g.coefficients()] or [1])
    h = g.parent()(den * g)
    content = gcd(*[abs(ZZ(c)) for c in h.coefficients()])
    return h if content <= 1 else h // content


def fiber_polynomials(ring, eps, r, b):
    M, U, v = ring.gens()
    a = (eps * (r - 1) ** 3 - 1 - r**3 - b * r**2) / r
    c1 = 2 * a + 6
    c2 = a**2 + 2 * b - 15
    c3 = 2 * a * b + 22
    c4 = 2 * a + b**2 - 15
    c5 = 2 * b + 6
    B = c5 * M + 3 * U
    delta = 4 * c4 * M + 12 * (U**2 + v**2) - B**2
    f3 = B * delta + 16 * v**3 - 8 * c3 * M - 8 * U**3 - 48 * U * v**2
    f2 = delta**2 + 64 * B * v**3 - 64 * c2 * M - 192 * (U**2 * v**2 + v**4)
    f1 = delta * v**3 - 4 * c1 * M - 12 * U * v**4
    return [primitive_integral(g) for g in (f1, f2, f3)], a


def auto_a(b, r):
    return [
        2*b**3*r**4 + 30*b**2*r**4 - 36*b**2*r**3 + 8*b**2*r**2
        + 126*b*r**4 - 296*b*r**3 + 216*b*r**2 - 48*b*r
        + 162*r**4 - 612*r**3 + 816*r**2 - 464*r + 96,
        -b**4*r**4 - 12*b**3*r**4 + 12*b**3*r**3 - 2*b**3*r**2
        - 54*b**2*r**4 + 108*b**2*r**3 - 66*b**2*r**2 + 12*b**2*r
        - 108*b*r**4 + 324*b*r**3 - 342*b*r**2 + 152*b*r - 24*b
        - 81*r**4 + 324*r**3 - 470*r**2 + 300*r - 72,
        b**4*r**3 + 6*b**3*r**3 - 6*b**3*r**2 - 14*b**2*r**2
        + 12*b**2*r - 54*b*r**3 + 102*b*r**2 - 48*b*r
        - 81*r**3 + 270*r**2 - 284*r + 96,
    ]


def auto_b(b, r):
    return [
        12*b**3*r**4 - 12*b**3*r**3 + 4*b**3*r**2 + 108*b**2*r**4
        - 188*b**2*r**3 + 120*b**2*r**2 - 36*b**2*r + 4*b**2
        + 324*b*r**4 - 804*b*r**3 + 732*b*r**2 - 296*b*r + 48*b
        + 324*r**4 - 1044*r**3 + 1224*r**2 - 612*r + 108,
        -6*b**3*r**4 + 6*b**3*r**3 - 2*b**3*r**2 - 54*b**2*r**4
        + 94*b**2*r**3 - 66*b**2*r**2 + 24*b**2*r - 4*b**2
        - 162*b*r**4 + 378*b*r**3 - 342*b*r**2 + 144*b*r - 24*b
        - 162*r**4 + 450*r**3 - 470*r**2 + 216*r - 36,
        6*b**2*r**2 - 6*b**2*r + 2*b**2 + 24*b*r**3 - 24*b*r**2
        + 4*b*r + 72*r**3 - 142*r**2 + 90*r - 18,
    ]


def exact_auto_label(b, r):
    if all(e == 0 for e in auto_a(b, r)):
        return "A"
    if all(e == 0 for e in auto_b(b, r)):
        return "B"
    return "none"


def good_auxiliary_root(eps, r, b, z, p):
    field = GF(p)
    rr = field(r.numerator()) / field(r.denominator())
    bb = field(b.numerator()) / field(b.denominator())
    ee = field(eps)
    mm, uu, vv = map(field, z)
    if not mm or not mm.is_square() or not vv or uu**2 - 4*vv**2 == 0:
        return False
    poly = PolynomialRing(field, "x")
    x = poly.gen()
    aa = (ee * (rr - 1) ** 3 - 1 - rr**3 - bb * rr**2) / rr
    h = 1 + aa*x + bb*x**2 + x**3
    f = h**2 - (x - 1) ** 6
    q = x**2 + uu*x + vv**2
    return (
        f.degree() == 5
        and f.discriminant() != 0
        and h(field(1)) != 0
        and f(rr) == 0
        and q.discriminant() != 0
        and q.gcd(f).degree() == 0
    )


def smooth_roots(polys, eps, r, b, p):
    field = GF(p)
    derivs = [[g.derivative(g.parent().gen(j)) for j in range(3)] for g in polys]
    roots = []
    for mm, uu, vv in product(range(1, p), range(p), range(1, p)):
        z = (mm, uu, vv)
        if any(ZZ(g(*z)) % p for g in polys):
            continue
        if not good_auxiliary_root(eps, r, b, z, p):
            continue
        jac = matrix(field, [[derivs[i][j](*z) for j in range(3)] for i in range(3)])
        if jac.det() != 0:
            roots.append(z)
    return roots


def hensel_lift(polys, z0, p, precision):
    derivs = [[g.derivative(g.parent().gen(j)) for j in range(3)] for g in polys]
    field = GF(p)
    z = [ZZ(c) for c in z0]
    modulus = ZZ(p)
    for _ in range(2, precision + 1):
        values = [ZZ(g(*z)) for g in polys]
        if any(value % modulus for value in values):
            return None
        jac = matrix(field, [[derivs[i][j](*z) for j in range(3)] for i in range(3)])
        if jac.det() == 0:
            return None
        rhs = vector(field, [-(value // modulus) for value in values])
        delta = jac.solve_right(rhs)
        new_modulus = modulus * p
        z = [ZZ(z[j] + modulus * ZZ(delta[j])) % new_modulus for j in range(3)]
        modulus = new_modulus
    if any(ZZ(g(*z)) % modulus for g in polys):
        return None
    return tuple(z)


def balanced_rational_reconstruction(a, modulus):
    aa = ZZ(a) % modulus
    if aa == 0:
        return QQ(0)
    bound = isqrt(ZZ(modulus // 2))
    r0, r1 = ZZ(modulus), aa
    t0, t1 = ZZ(0), ZZ(1)
    while abs(r1) > bound:
        quotient = r0 // r1
        r0, r1 = r1, r0 - quotient * r1
        t0, t1 = t1, t0 - quotient * t1
    if not t1 or abs(t1) > bound or abs(r1) > bound:
        return None
    num, den = r1, t1
    if den < 0:
        num, den = -num, -den
    if gcd(num, den) != 1 or (num - aa * den) % modulus:
        return None
    return QQ(num) / den


def exact_fiber_analysis(log, label, ring, polys, do_primary):
    M, U, v = ring.gens()
    ideal = ring.ideal(polys)
    saturated, saturation_exponent = ideal.saturation(ring.ideal([M*v*(U**2 - 4*v**2)]))
    log(f"EXACT_FIBER {label}")
    log(
        f" saturated_dimension {saturated.dimension()} degree {saturated.degree()} "
        f"saturation_exponent {saturation_exponent}"
    )
    basis = saturated.groebner_basis()
    log(f" groebner_length {len(basis)}")
    log(f" degree_term_summary {[(g.total_degree(), len(g.monomials())) for g in basis]}")
    univariate_v = [g for g in basis if g.degree(M) == 0 and g.degree(U) == 0]
    log(f" univariate_v_count {len(univariate_v)}")
    for g in univariate_v:
        factors = [(factor.degree(v), exponent) for factor, exponent in g.factor()]
        log(f"  univariate_v_degree {g.degree(v)} factor_degrees {factors}")
        log(f"  univariate_v {g}")
    if do_primary:
        components = saturated.primary_decomposition()
        log(f" primary_components {len(components)}")
        for index, component in enumerate(components, 1):
            log(
                f"  component {index} dimension {component.dimension()} "
                f"degree {component.degree()} prime {component.is_prime()}"
            )
    return saturated


def try_crt(log, label, polys, p1, lifts1, p2, lifts2, precision):
    m1, m2 = ZZ(p1**precision), ZZ(p2**precision)
    modulus = m1 * m2
    pairs = reconstructed = exact = 0
    seen = set()
    for z1, z2 in product(lifts1, lifts2):
        pairs += 1
        values = []
        for index in range(3):
            residue = crt(ZZ(z1[index]), ZZ(z2[index]), m1, m2)
            value = balanced_rational_reconstruction(residue, modulus)
            if value is None:
                break
            values.append(value)
        if len(values) != 3:
            continue
        reconstructed += 1
        key = tuple(values)
        if key in seen:
            continue
        seen.add(key)
        if all(g(*values) == 0 for g in polys):
            exact += 1
            log(
                f" EXACT_CRT_POINT {label} {values} "
                f"M_square {values[0].is_square()}"
            )
    log(
        f"CRT_SUMMARY {label} primes {p1} {p2} precision {precision} "
        f"modulus {modulus} balanced_height_bound {isqrt(modulus // 2)} "
        f"pairs {pairs} reconstructed {reconstructed} exact_points {exact}"
    )


def main():
    parser = ArgumentParser()
    parser.add_argument("--lift-precision", type=int, default=4)
    parser.add_argument("--primary", action="store_true")
    parser.add_argument("--skip-exact", action="store_true")
    parser.add_argument(
        "--output", default="data/contact6_m36_266_targeted_lift.txt"
    )
    args = parser.parse_args()
    if args.lift_precision < 2 or args.lift_precision > 8:
        raise ValueError("lift precision must be between 2 and 8")

    lines = []

    def log(message):
        text = str(message)
        print(text, flush=True)
        lines.append(text)

    ring = PolynomialRing(QQ, names=("M", "U", "v"), order="lex")
    fibers = [
        ("r2_b3", QQ(1), QQ(2), QQ(3), [13, 19]),
        ("r4_b5", QQ(1), QQ(4), QQ(5), [11, 13, 17, 19, 23, 29, 31]),
    ]
    log("Bounded [2,6,6] component/Hensel/CRT diagnostic (Sage/Singular)")
    log(f"lift_precision {args.lift_precision} do_primary {args.primary} skip_exact {args.skip_exact}")
    log("go_stop all smooth roots; fixed primes <=31; no height box")

    for label, eps, r, b, primes in fibers:
        polys, a = fiber_polynomials(ring, eps, r, b)
        log(f"\nFIBER {label} eps {eps} r {r} a {a} b {b}")
        log(
            f" exact_auto_label {exact_auto_label(b, r)} "
            f"AutoA {auto_a(b, r)} AutoB {auto_b(b, r)}"
        )
        log(
            " polynomial_summaries "
            + str([(g.total_degree(), len(g.monomials())) for g in polys])
        )
        (log(" exact_fiber_analysis SKIPPED_BY_TIME_CAP") if args.skip_exact else exact_fiber_analysis(log, label, ring, polys, args.primary))

        records = {}
        for p in primes:
            if a.denominator() % p == 0:
                continue
            roots = smooth_roots(polys, eps, r, b, p)
            lifts = [
                lifted
                for root in roots
                if (lifted := hensel_lift(polys, root, p, args.lift_precision))
                is not None
            ]
            automorphic_mod_p = all(ZZ(e.numerator()) % p == 0 for e in auto_a(b, r)) or all(
                ZZ(e.numerator()) % p == 0 for e in auto_b(b, r)
            )
            log(
                f" PRIME {p} smooth_good_roots {len(roots)} lifted {len(lifts)} "
                f"auto_mod_p {automorphic_mod_p} root_samples {roots[:4]}"
            )
            records[p] = lifts

        if label == "r2_b3":
            try_crt(log, label, polys, 13, records[13], 19, records[19], args.lift_precision)
        else:
            second_primes = [p for p in primes if p != 19 and records.get(p)]
            if second_primes:
                p2 = second_primes[0]
                try_crt(log, label, polys, p2, records[p2], 19, records[19], args.lift_precision)
            else:
                log(f"CRT_SUMMARY {label} no second useful prime in bounded scan")

    log("\nDONE")
    with open(args.output, "w", encoding="ascii") as handle:
        handle.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
