#!/usr/bin/env python3
"""Symbolic diagnostics for the reduced [8,8] square subcover.

This derives the small eliminants used after the height-80 reduced search:

* the first-cover equation on the square branch, factored into two U-branches;
* the universal tau-elimination for the second halving;
* the M,N sign-eliminated second-halving condition H_eta(R,s,U,z).

The script intentionally prints summaries instead of full 10k-term
polynomials.  The full U-resultant against a first-cover branch is much larger
than the useful intermediate equations and is not materialized here.
"""

import sympy as sp


def summary(label, poly, gens):
    poly = sp.Poly(poly, *gens)
    print(
        label,
        "total_degree",
        poly.total_degree(),
        "degree_by_var",
        [poly.degree(g) for g in gens],
        "terms",
        len(poly.terms()),
    )


def main():
    U, M, N, c, z, mu, c1, c2, c3 = sp.symbols("U M N c z mu c1 c2 c3")

    a = 2 * mu - U
    d = a - c * z**2
    pz = (
        c**2 * z**4
        - 4 * M * c * z**3
        + (4 * M**2 + 2 * U * c) * z**2
        + (4 * M * U - 8 * N) * z
        + U**2
        - 4 * mu**2
    )
    b = 2 * z * (M * mu - N) - mu * a - c * mu * z**2

    print("Reduced second-halving equations after solving the x^3 coefficient:")
    print("  E2 = (P(z) + 4*D(z)*tau)/4")
    print("  E1 = (2*mu-U)*tau^2 + B(z)*tau + mu^2*c*z^2")
    print("  D(z) =", d)
    print("  P(z) =", pz)
    print("  B(z) =", b)

    # Eliminating tau without dividing by D gives:
    #   G = A*P^2 - 4*B*P*D + 16*mu^2*c*z^2*D^2.
    g = sp.expand(a * pz**2 - 4 * b * pz * d + 16 * mu**2 * c * z**2 * d**2)
    summary("UNIVERSAL_G", g, [U, M, N, c, mu, z])

    # First-cover relations:
    #   M^2 = c3 - 2*c*U,
    #   N^2 = c1 - 2*c*U*mu^2,
    #   2MN = c2 - c*(U^2 + 2*mu^2).
    m2 = c3 - 2 * c * U
    n2 = c1 - 2 * c * U * mu**2
    mn2 = c2 - c * (U**2 + 2 * mu**2)

    gb = sp.groebner(
        [M**2 - m2, N**2 - n2, 2 * M * N - mn2],
        M,
        N,
        z,
        U,
        c,
        mu,
        c1,
        c2,
        c3,
        order="lex",
    )
    rem = sp.expand(gb.reduce(g)[1])
    lm = sp.expand(rem.coeff(M))
    ln = sp.expand(rem.coeff(N))
    l0 = sp.expand(rem.subs({M: 0, N: 0}))
    print()
    print("Modulo the first-cover relations, G is linear in M,N:")
    summary("  L_M", lm, [U, c, mu, c1, c2, c3, z])
    summary("  L_N", ln, [U, c, mu, c1, c2, c3, z])
    summary("  L_0", l0, [U, c, mu, c1, c2, c3, z])

    h_generic = sp.expand(l0**2 - lm**2 * m2 - ln**2 * n2 - lm * ln * mn2)
    summary("GENERIC_SIGN_ELIMINANT_H", h_generic, [U, c, mu, c1, c2, c3, z])

    R, s, x = sp.symbols("R s x")
    w = s**2
    t = (2 * R**2 + (1 - w**2) * R - 2 * w**2) / (4 * (w**2 - 1))
    a_poly = x**2 + (R**3 + 4 * R**2 * t + R - 8 * R * t + 4 * t) * x + R**4
    b_poly = (
        (R + 2 + 4 * t) * x**2
        + (R**2 + 4 * R + 1 + 8 * t) * x
        + (2 * R**2 + R + 4 * t)
    )
    h = sp.expand(a_poly * b_poly)
    coeffs = {
        c1: sp.factor(h.coeff(x, 1)),
        c2: sp.factor(h.coeff(x, 2)),
        c3: sp.factor(h.coeff(x, 3)),
        c: sp.factor(h.coeff(x, 4)),
    }

    first_m2 = coeffs[c3] - 2 * coeffs[c] * U
    first_n2 = coeffs[c1] - 2 * coeffs[c] * U * (R * s) ** 2
    first_mn2 = coeffs[c2] - coeffs[c] * (U**2 + 2 * (R * s) ** 2)
    first = sp.factor(sp.together(4 * first_m2 * first_n2 - first_mn2**2).as_numer_denom()[0])

    print()
    print("First-cover equation on w=s^2:")
    for i, (fac, exp) in enumerate(sp.factor_list(first)[1], 1):
        poly = sp.Poly(fac, R, s, U)
        print(
            "  factor",
            i,
            "exp",
            exp,
            "total_degree",
            poly.total_degree(),
            "degree_U",
            poly.degree(U),
            "terms",
            len(poly.terms()),
        )

    print()
    print("Substituted sign-eliminated second-cover equations:")
    for sign in [1, -1]:
        h_eta = sp.together(h_generic.subs(coeffs).subs(mu, sign * R * s)).as_numer_denom()[0]
        h_eta = sp.factor(h_eta)
        summary(f"  H_eta_{sign:+d}", h_eta, [R, s, U, z])
        factors = sp.factor_list(h_eta)[1]
        print("    irreducible_factor_count", len(factors))
        for fac, exp in factors:
            poly = sp.Poly(fac, R, s, U, z)
            print(
                "    factor_exp",
                exp,
                "total_degree",
                poly.total_degree(),
                "degree_U",
                poly.degree(U),
                "degree_z",
                poly.degree(z),
                "terms",
                len(poly.terms()),
            )

    print()
    print("U-resultant note:")
    print("  Each nonboundary first-cover branch has degree_U 2.")
    print("  H_eta has degree_U 10 and degree_z 16 after substitution.")
    print("  A direct resultant Res_U(branch,H_eta) was tested and is too large")
    print("  for routine interactive use; finite-prime diagnostics should use")
    print("  the reduced E1/E2 formulas instead.")


if __name__ == "__main__":
    main()
