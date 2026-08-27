#!/usr/bin/env sage
"""Bounded component diagnostics for the non-automorphic [2,6,6] surface.

The eps=+1 extra-root chart has coordinates (r,b,M,U,v), with M=L^2.
Three cubic-contact equations cut out an expected surface.  This script
compares the raw ideal, its open-chart saturation, and saturation away from
the two explicit extra-involution ideals.  It also counts smooth points on
the corresponding open set over small finite fields.
"""

from argparse import ArgumentParser
from sage.all import GF, QQ, PolynomialRing, factorial, matrix


def automorphism_polynomials(b, r):
    auto_a = [
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
    auto_b = [
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
    return auto_a, auto_b


def chart(field):
    ring = PolynomialRing(field, names=("r", "b", "M", "U", "v"),
                          order="degrevlex")
    r, b, M, U, v = ring.gens()
    fraction = ring.fraction_field()
    x = PolynomialRing(fraction, "x").gen()
    a = 3 - (b + 3)*r - 2/r
    h = 1 + a*x + b*x**2 + x**3
    f = h**2 - (x - 1)**6
    coeffs = [f[i] for i in range(6)]
    B = coeffs[5]*M + 3*U
    delta = 4*coeffs[4]*M + 12*(U**2 + v**2) - B**2
    equations = [
        delta*v**3 - 4*coeffs[1]*M - 12*U*v**4,
        delta**2 + 64*B*v**3 - 64*coeffs[2]*M
        - 192*(U**2*v**2 + v**4),
        B*delta + 16*v**3 - 8*coeffs[3]*M - 8*U**3 - 48*U*v**2,
    ]
    equations = [ring(g.numerator()) for g in equations]
    q = x**2 + U*x + v**2
    smooth_f = ring(f.discriminant().numerator())
    coprime_qf = ring(f.resultant(q).numerator())
    h_at_one = ring((a + b + 2).numerator())
    structural = r*(r - 1)*(b + 3)*M*v*(U**2 - 4*v**2)*h_at_one
    auto_a, auto_b = automorphism_polynomials(b, r)
    return {
        "ring": ring,
        "gens": (r, b, M, U, v),
        "equations": equations,
        "structural": structural,
        "smooth_f": smooth_f,
        "coprime_qf": coprime_qf,
        "auto_a": auto_a,
        "auto_b": auto_b,
        "g_auto_a": auto_a[0].gcd(auto_a[1]).gcd(auto_a[2]),
        "g_auto_b": auto_b[0].gcd(auto_b[1]).gcd(auto_b[2]),
    }


def fixed_r_chart(data, fixed_r):
    old_ring = data["ring"]
    old_r, old_b, old_M, old_U, old_v = old_ring.gens()
    field = old_ring.base_ring()
    ring = PolynomialRing(field, names=("b", "M", "U", "v"),
                          order="degrevlex")
    b, M, U, v = ring.gens()
    phi = old_ring.hom([field(fixed_r), b, M, U, v], ring)
    return {
        "ring": ring,
        "gens": (b, M, U, v),
        "equations": [phi(g) for g in data["equations"]],
        "structural": phi(data["structural"]),
        "smooth_f": phi(data["smooth_f"]),
        "coprime_qf": phi(data["coprime_qf"]),
        "auto_a": [phi(g) for g in data["auto_a"]],
        "auto_b": [phi(g) for g in data["auto_b"]],
        "g_auto_a": phi(data["g_auto_a"]),
        "g_auto_b": phi(data["g_auto_b"]),
    }


def ideal_summary(log, label, ideal, with_degree=True):
    log(f"IDEAL {label} dimension {ideal.dimension()}")
    if with_degree:
        try:
            log(f"IDEAL {label} degree {ideal.degree()}")
        except Exception as error:
            log(f"IDEAL {label} degree ERROR {type(error).__name__}: {error}")
    basis = ideal.groebner_basis()
    log(f"IDEAL {label} groebner_length {len(basis)}")
    log("IDEAL %s groebner_summary %s" % (
        label,
        [(g.total_degree(), len(g.monomials())) for g in basis],
    ))
    return basis


def short_ideal_summary(log, label, ideal):
    log(f"IDEAL {label} dimension {ideal.dimension()}")
    basis = ideal.groebner_basis()
    log(f"IDEAL {label} groebner_length {len(basis)}")
    log("IDEAL %s groebner_summary %s" % (
        label,
        [(g.total_degree(), len(g.monomials())) for g in basis],
    ))
    try:
        dimension = ideal.dimension()
        projective_closure = ideal.homogenize()
        hilbert = projective_closure.hilbert_polynomial()
        degree = hilbert.leading_coefficient()*factorial(dimension)
        log(f"IDEAL {label} hilbert_polynomial {hilbert} affine_degree {degree}")
    except Exception as error:
        log(f"IDEAL {label} hilbert_degree ERROR {type(error).__name__}: {error}")
    return basis


def saturated(log, label, ideal, away):
    result, exponent = ideal.saturation(away)
    log(f"SATURATION {label} exponent {exponent}")
    return result


def algebraic_diagnostic(log, data, decompose, full_open):
    ring = data["ring"]
    ideal = ring.ideal(data["equations"])
    log(f"equation_summaries {[(g.total_degree(), len(g.monomials())) for g in data['equations']]}")
    log(f"structural_boundary degree {data['structural'].total_degree()} terms {len(data['structural'].monomials())}")
    log(f"smooth_discriminant degree {data['smooth_f'].total_degree()} terms {len(data['smooth_f'].monomials())}")
    log(f"q_f_resultant degree {data['coprime_qf'].total_degree()} terms {len(data['coprime_qf'].monomials())}")
    log(f"auto_gcd_A {data['g_auto_a']}")
    log(f"auto_gcd_B {data['g_auto_b']}")
    raw_basis = ideal_summary(log, "raw", ideal)
    structural = saturated(log, "structural", ideal, ring.ideal([data["structural"]]))
    structural_basis = ideal_summary(log, "structural", structural)
    log(f"COMPARE raw_equals_structural {ring.ideal(raw_basis) == ring.ideal(structural_basis)}")
    opened = structural
    if full_open:
        opened = saturated(log, "smooth_curve", opened, ring.ideal([data["smooth_f"]]))
        opened = saturated(log, "coprime_q_f", opened, ring.ideal([data["coprime_qf"]]))
        ideal_summary(log, "full_open", opened)
    no_a = saturated(log, "away_auto_A", opened, ring.ideal(data["auto_a"]))
    no_ab = saturated(log, "away_auto_B", no_a, ring.ideal(data["auto_b"]))
    no_ab_basis = ideal_summary(log, "nonautomorphic_open", no_ab)
    log(f"COMPARE open_equals_nonautomorphic {opened == ring.ideal(no_ab_basis)}")
    # Principal gcds define the divisorial parts of the two automorphism loci.
    for label, g in (("A", data["g_auto_a"]), ("B", data["g_auto_b"])):
        auto_section = opened + ring.ideal([g])
        ideal_summary(log, f"auto_divisor_{label}", auto_section)
    if decompose:
        components = no_ab.primary_decomposition()
        log(f"PRIMARY count {len(components)}")
        for index, component in enumerate(components, 1):
            log(
                f"PRIMARY {index} dimension {component.dimension()} "
                f"degree {component.degree()} prime {component.is_prime()}"
            )


def fixed_r_diagnostic(log, data, decompose, full_open, prime_test):
    ring = data["ring"]
    ideal = ring.ideal(data["equations"])
    log(f"fixed_equation_summaries {[(g.total_degree(), len(g.monomials())) for g in data['equations']]}")
    raw_basis = short_ideal_summary(log, "fixed_raw", ideal)
    structural = saturated(log, "fixed_structural", ideal,
                           ring.ideal([data["structural"]]))
    structural_basis = short_ideal_summary(log, "fixed_structural", structural)
    log(f"COMPARE fixed_raw_equals_structural {ring.ideal(raw_basis) == ring.ideal(structural_basis)}")
    opened = structural
    if full_open:
        opened = saturated(log, "fixed_smooth_curve", opened,
                           ring.ideal([data["smooth_f"]]))
        opened = saturated(log, "fixed_coprime_q_f", opened,
                           ring.ideal([data["coprime_qf"]]))
        short_ideal_summary(log, "fixed_full_open", opened)
    no_a = saturated(log, "fixed_away_auto_A", opened,
                     ring.ideal(data["auto_a"]))
    no_ab = saturated(log, "fixed_away_auto_B", no_a,
                      ring.ideal(data["auto_b"]))
    no_ab_basis = short_ideal_summary(log, "fixed_nonautomorphic_open", no_ab)
    log(f"COMPARE fixed_open_equals_nonautomorphic {opened == ring.ideal(no_ab_basis)}")
    if prime_test:
        log(f"FIXED_IS_PRIME {no_ab.is_prime()}")
    if decompose:
        components = no_ab.primary_decomposition()
        log(f"FIXED_PRIMARY count {len(components)}")
        for index, component in enumerate(components, 1):
            basis = component.groebner_basis()
            log(
                f"FIXED_PRIMARY {index} dimension {component.dimension()} "
                f"prime {component.is_prime()} groebner_summary "
                f"{[(g.total_degree(), len(g.monomials())) for g in basis]}"
            )


def evaluate_good_point(data, values):
    ring = data["ring"]
    point = dict(zip(ring.gens(), values))
    if any(g.subs(point) != 0 for g in data["equations"]):
        return None
    if data["structural"].subs(point) == 0:
        return None
    if data["smooth_f"].subs(point) == 0 or data["coprime_qf"].subs(point) == 0:
        return None
    auto_a = all(g.subs(point) == 0 for g in data["auto_a"])
    auto_b = all(g.subs(point) == 0 for g in data["auto_b"])
    return auto_a, auto_b


def finite_count(log, data):
    ring = data["ring"]
    field = ring.base_ring()
    p = field.cardinality()
    equations = data["equations"]
    jacobian = [[g.derivative(x) for x in ring.gens()] for g in equations]
    total = good = nonauto = auto_a = auto_b = both = smooth_nonauto = 0
    square_m_nonauto = 0
    base_nonauto = set()
    samples = []
    # The open conditions already require r,M,v nonzero and r != 1.
    for rr in field:
        if rr == 0 or rr == 1:
            continue
        for bb in field:
            for mm in field:
                if mm == 0:
                    continue
                for uu in field:
                    for vv in field:
                        if vv == 0:
                            continue
                        values = (rr, bb, mm, uu, vv)
                        if any(g(*values) != 0 for g in equations):
                            continue
                        total += 1
                        kind = evaluate_good_point(data, values)
                        if kind is None:
                            continue
                        good += 1
                        in_a, in_b = kind
                        auto_a += int(in_a)
                        auto_b += int(in_b)
                        both += int(in_a and in_b)
                        if in_a or in_b:
                            continue
                        nonauto += 1
                        square_m_nonauto += int(mm.is_square())
                        base_nonauto.add((int(rr), int(bb)))
                        jac = matrix(field, [[entry(*values) for entry in row] for row in jacobian])
                        tangent_dimension = 5 - jac.rank()
                        smooth_nonauto += int(tangent_dimension == 2)
                        if len(samples) < 8:
                            samples.append((tuple(map(int, values)), tangent_dimension))
    log(
        f"FINITE p {p} equation_points_open_loop {total} good {good} "
        f"auto_A {auto_a} auto_B {auto_b} both {both} nonauto {nonauto}"
    )
    log(
        f"FINITE p {p} smooth_nonauto {smooth_nonauto} "
        f"square_M_nonauto {square_m_nonauto} distinct_nonauto_bases {len(base_nonauto)}"
    )
    log(f"FINITE p {p} nonauto_samples {samples}")


def main():
    parser = ArgumentParser()
    parser.add_argument("--prime", type=int, default=13)
    parser.add_argument("--decompose", action="store_true")
    parser.add_argument("--full-open", action="store_true")
    parser.add_argument("--count", action="store_true")
    parser.add_argument("--count-only", action="store_true")
    parser.add_argument("--fixed-r", type=int, default=-1)
    parser.add_argument("--rational", action="store_true")
    parser.add_argument("--prime-test", action="store_true")
    parser.add_argument("--output", default="")
    args = parser.parse_args()
    lines = []

    def log(message):
        message = str(message)
        print(message, flush=True)
        lines.append(message)

    log("Contact-6 eps=+1 [2,6,6] surface decomposition diagnostic")
    log(
        f"prime {args.prime} decompose {args.decompose} full_open {args.full_open} "
        f"count {args.count} count_only {args.count_only} fixed_r {args.fixed_r} "
        f"rational {args.rational}"
    )
    data = chart(QQ if args.rational else GF(args.prime))
    if not args.count_only:
        if args.fixed_r >= 0:
            fixed = fixed_r_chart(data, args.fixed_r)
            fixed_r_diagnostic(log, fixed, args.decompose, args.full_open,
                               args.prime_test)
        else:
            algebraic_diagnostic(log, data, args.decompose, args.full_open)
    if args.count and not args.rational:
        finite_count(log, data)
    output = args.output or f"data/contact6_m36_266_surface_decompose_p{args.prime}.txt"
    with open(output, "w", encoding="ascii") as handle:
        handle.write("\n".join(lines) + "\n")
    log(f"wrote {output}")


if __name__ == "__main__":
    main()
