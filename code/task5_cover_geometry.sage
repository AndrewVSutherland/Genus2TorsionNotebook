#!/usr/bin/env sage
"""Bounded geometry triage for the reduced [8,8] and corrected [64] covers.

This script performs exact symbolic computations only.  It deliberately does
not enumerate rational parameters.  Expensive invocations should be wrapped in
an external ``timeout``; the checkpoints make partial progress reproducible.
"""

from sage.all import *


def poly_summary(label, f, variables):
    """Print stable degree/term data for a polynomial or rational numerator."""
    if f in (0, 1, -1):
        print(f"{label}: {f}")
        return
    p = f if hasattr(f, "dict") else f.numerator()
    parent_gens = {str(g): g for g in p.parent().gens()}
    summary_gens = [parent_gens[str(v)] for v in variables]
    print(
        f"{label}: total_degree={p.total_degree()} "
        f"degrees={[p.degree(v) for v in summary_gens]} terms={len(p.dict())}"
    )


def c32_polynomial(z, r, w):
    """Homogeneous degree-8 closure of the reconstructed Elkies [32] base."""
    return (
        3*z**4*r**4 + 9*z**4*r**3*w - 16*z**3*r**4*w
        + 10*z**4*r**2*w**2 - 56*z**3*r**3*w**2 + 32*z**2*r**4*w**2
        + 5*z**4*r*w**3 - 72*z**3*r**2*w**3 + 144*z**2*r**3*w**3
        - 64*z*r**4*w**3 + z**4*w**4 - 40*z**3*r*w**4
        + 208*z**2*r**2*w**4 - 224*z*r**3*w**4 + 80*r**4*w**4
        - 8*z**3*w**5 + 120*z**2*r*w**5 - 288*z*r**2*w**5
        + 160*r**3*w**5 + 24*z**2*w**6 - 160*z*r*w**6
        + 160*r**2*w**6 - 32*z*w**7 + 80*r*w**7 + 16*w**8
    )


def p_num(z, r):
    return (
        (z - 2)*(r + 1)**2*(z**2*r + z**2 - 4*r)
        * (z**4*r**2 + 2*z**4*r - 10*z**3*r**2 + z**4
           - 16*z**3*r + 28*z**2*r**2 - 6*z**3 + 44*z**2*r
           - 16*z*r**2 + 16*z**2 - 56*z*r + 16*r**2
           - 24*z + 32*r + 16)
    )


def p_den(z, r):
    return r*z*(
        z**8*r**4 + 4*z**8*r**3 - 24*z**7*r**4 + 48*z**6*r**5
        + 6*z**8*r**2 - 80*z**7*r**3 + 284*z**6*r**4 - 256*z**5*r**5
        + 4*z**8*r - 100*z**7*r**2 + 632*z**6*r**3 - 1200*z**5*r**4
        + 320*z**4*r**5 + z**8 - 56*z**7*r + 672*z**6*r**2
        - 2432*z**5*r**3 + 2224*z**4*r**4 - 12*z**7 + 344*z**6*r
        - 2496*z**5*r**2 + 5312*z**4*r**3 - 2304*z**3*r**4 + 68*z**6
        - 1248*z**5*r + 5808*z**4*r**2 - 7232*z**3*r**3 + 1664*z**2*r**4
        - 240*z**5 + 2976*z**4*r - 8832*z**3*r**2 + 6272*z**2*r**3
        - 768*z*r**4 + 576*z**4 - 4800*z**3*r + 8640*z**2*r**2
        - 3328*z*r**3 + 256*r**4 - 960*z**3 + 5120*z**2*r
        - 5120*z*r**2 + 1024*r**3 + 1088*z**2 - 3328*z*r
        + 1536*r**2 - 768*z + 1024*r + 256
    )


def small_c32_parameter(t):
    """Clean normalization coordinate obtained from the factored degree-8 map."""
    r = 1/(t**4 - 1)
    z = 2*(t**4 + t**3 + t**2 + t + 1)/(t**2*(t**2 + t + 1))
    return z, r


def analyze_c32_base():
    print("=== corrected [64]: normalize/parametrize the [32] base ===")
    P2 = ProjectiveSpace(QQ, 2, names=("Z", "R", "W"))
    Z, R, W = P2.gens()
    C = Curve(c32_polynomial(Z, R, W))
    print(f"base_projective_degree={C.degree()}")
    print(f"base_arithmetic_genus={C.arithmetic_genus()}")
    print(f"base_geometric_genus={C.geometric_genus()}")
    print(f"base_irreducible={C.is_irreducible()}")
    print("computing_rational_parameterization")
    phi = C.rational_parameterization()
    print(f"parameterization={phi}")
    try:
        eqs = phi.defining_polynomials()
    except (AttributeError, TypeError):
        eqs = phi
    print(f"parameterization_coordinate_degrees={[e.degree() for e in eqs]}")
    for i, e in enumerate(eqs):
        print(f"parameterization_coordinate_{i}={e}")


def analyze_halving_generic():
    """Eliminate one halving variable for a generic normalized quintic."""
    print("=== generic point-halving fiber ===")
    R = PolynomialRing(QQ, names=("f1", "f2", "f3", "f4", "f5", "m", "n"))
    f1, f2, f3, f4, f5, m, n = R.gens()
    n_polynomial = n
    K = R.fraction_field()
    f1, f2, f3, f4, f5, m, n = map(K, (f1, f2, f3, f4, f5, m, n))
    # k=+1 suffices: k=-1 is carried to it by (m,n)->(-m,-n).
    alpha = (f4 - m**2)/(2*f5)
    beta = (f3 - 2*m*n - f5*alpha**2)/(2*f5)
    k2 = f2 - n**2 - 2*m - 2*f5*alpha*beta
    k1 = f1 - 2*n - f5*beta**2
    k2n = R(k2.numerator())
    k1n = R(k1.numerator())
    poly_summary("K2_generic", k2n, R.gens())
    poly_summary("K1_generic", k1n, R.gens())
    print("computing_generic_resultant_in_n")
    res = k1n.resultant(k2n, n_polynomial)
    fac = res.factor()
    poly_summary("Res_n(K1,K2)", res, (f1, f2, f3, f4, f5, m))
    print(f"resultant_factor_count={len(fac)}")
    for i, (g, exponent) in enumerate(fac):
        print(f"resultant_factor_{i}_exponent={exponent}")
        poly_summary(f"resultant_factor_{i}", g, (f1, f2, f3, f4, f5, m))


def analyze_64_extension():
    """Build the corrected order-64 function field over the clean P1 base."""
    print("=== corrected [64] degree-16 function field ===")
    K = FunctionField(QQ, "t")
    t = K.gen()
    z, r = small_c32_parameter(t)
    p = 4*p_num(z, r)/p_den(z, r)
    a = z*p
    c = (-p*z**4 + 4*p*z**3 - 8*p*z**2 + 8*z**2 + 8*p*z - 16)/(4*z**2)
    b = p*(z - 1) + c - 1
    print(f"z={z.factor()}")
    print(f"r={r.factor()}")
    print(f"p={p.factor()}")
    print(f"a={a.factor()}")
    print(f"b={b.factor()}")
    print(f"c={c.factor()}")
    fvalues = [2*c, c**2 + 2*b, 2*a + 2*b*c, b**2 + 2*a*c, a*(2*b - a)]
    PR = PolynomialRing(QQ, names=("f1", "f2", "f3", "f4", "f5", "m", "n"))
    f1, f2, f3, f4, f5, m, n = PR.gens()
    PK = PR.fraction_field()
    F1, F2, F3, F4, F5, M, N = map(PK, PR.gens())
    alpha = (F4 - M**2)/(2*F5)
    beta = (F3 - 2*M*N - F5*alpha**2)/(2*F5)
    k2 = PR((F2 - N**2 - 2*M - 2*F5*alpha*beta).numerator())
    k1 = PR((F1 - 2*N - F5*beta**2).numerator())
    H = [g for g, e in k1.resultant(k2, n).factor() if g.degree(m) == 16][0]
    KM = PolynomialRing(K, "M")
    mm = KM.gen()
    h = PR.hom(fvalues + [mm, K(0)], KM)(H).monic()
    print(f"extension_polynomial_degree={h.degree()}")
    print(f"extension_factor_degrees_exponents={[(g.degree(), e) for g, e in h.factor()]}")
    print(f"coefficient_max_numerator_degree={max(q.numerator().degree() for q in h.list())}")
    print(f"coefficient_max_denominator_degree={max(q.denominator().degree() for q in h.list())}")
    print("normalization_genus=not_computed (Sage FunctionField.genus timed out at 180s)")


def analyze_88_generic_degree():
    """Record the finite degrees visible before any large base elimination."""
    print("=== reduced [8,8] generic second-halving fiber ===")
    R = PolynomialRing(QQ, names=("U", "M", "N", "c", "z", "mu", "tau"))
    U, M, N, c, z, mu, tau = R.gens()
    A = 2*mu - U
    D = A - c*z**2
    P = (
        c**2*z**4 - 4*M*c*z**3 + (4*M**2 + 2*U*c)*z**2
        + (4*M*U - 8*N)*z + U**2 - 4*mu**2
    )
    B = 2*z*(M*mu - N) - mu*A - c*mu*z**2
    E2 = P + 4*D*tau
    E1 = A*tau**2 + B*tau + mu**2*c*z**2
    G = E1.resultant(E2, tau)
    poly_summary("E2", E2, (U, M, N, c, mu, z, tau))
    poly_summary("E1", E1, (U, M, N, c, mu, z, tau))
    poly_summary("Res_tau(E1,E2)", G, (U, M, N, c, mu, z))
    print(f"generic_degree_in_z={G.degree(z)}")
    print(f"generic_tau_degrees=({E1.degree(tau)},{E2.degree(tau)})")
    print("eta_branches=2; each has generically 8 z-roots and unique tau")


if __name__ == "__main__":
    import sys

    modes = sys.argv[1:] or ["base", "generic", "88"]
    if "base" in modes:
        analyze_c32_base()
    if "64" in modes:
        analyze_64_extension()
    if "generic" in modes:
        analyze_halving_generic()
    if "88" in modes:
        analyze_88_generic_degree()
