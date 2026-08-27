#!/usr/bin/env sage -python
"""
Reconstruct the algebraic conditions behind Elkies' N=32 genus-2 family.

The input from Elkies' page is one printed curve.  We normalize it into the
contact model

    y^2 = h(x)^2 - a^2*x^5*(x+1),      h = a*x^3 + b*x^2 + c*x + 1.

Then P0=(0,1), P1=(-1,h(-1)), and div(y-h)=5*P0+P1-6*oo, so
P1-oo = -5*(P0-oo).  If W=(r,0) is 16*(P0-oo), then

    4*P0 + 4*P1 + W - 9*oo

is principal.  Since L(9*oo) is spanned by
1,x,x^2,x^3,x^4,y,xy,x^2y, this is an explicit 9-by-8 rank drop.

This script builds those equations over Q, eliminates the linear contact
parameter, and records the component containing Elkies' printed example.
"""

from __future__ import annotations

from pathlib import Path

from sage.all import Matrix, PolynomialRing, QQ, gcd, prod


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "elkies32_reconstruct_conditions.txt"


def short_factor_line(fac, max_chars=1600):
    parts = []
    for base, exp in fac:
        s = str(base)
        if len(s) > max_chars:
            s = s[:max_chars] + " ...[truncated]"
        parts.append(f"({s})^{exp}")
    return " * ".join(parts) if parts else "1"


def coeff_by_var(poly, var_index, exp, target_ring):
    """Coefficient of the variable at var_index with exponent exp."""
    out = target_ring(0)
    for mon, coeff in poly.dict().items():
        if mon[var_index] == exp:
            new_mon = tuple(e for j, e in enumerate(mon) if j != var_index)
            out += target_ring({new_mon: coeff})
    return out


def lcm_poly(u, v):
    if u == 0:
        return v
    if v == 0:
        return u
    return (u * v) // gcd(u, v)


def main():
    lines = []
    lines.append("# Elkies N=32 reconstruction conditions")
    lines.append("")

    R = PolynomialRing(QQ, ("a", "b", "c", "r"))
    a, b, c, r = R.gens()
    K = R.fraction_field()

    X = PolynomialRing(K, "x")
    x = X.gen()

    h = K(a) * x**3 + K(b) * x**2 + K(c) * x + 1
    f = h**2 - K(a) ** 2 * x**5 * (x + 1)

    # L(9*oo) basis coefficients:
    # A(x)=A0+...+A4*x^4, B(x)=B0+B1*x+B2*x^2.
    # Rows are linear equations in [A0,A1,A2,A3,A4,B0,B1,B2].
    rows = []

    # Vanishing to order 4 at P0=(0,1).  The branch agrees with h to this
    # order because f-h^2 has an x^5 factor.
    h_coeffs_0 = [K(1), K(c), K(b), K(a)]
    for k in range(4):
        row = [K(0)] * 8
        if k <= 4:
            row[k] += 1
        for j in range(3):
            idx = k - j
            if 0 <= idx < len(h_coeffs_0):
                row[5 + j] += h_coeffs_0[idx]
        rows.append(row)

    # Vanishing to order 4 at P1=(-1,h(-1)).  Here we need the actual
    # Taylor expansion of the branch y=sqrt(f) through degree 3.
    T = PolynomialRing(K, "t")
    t = T.gen()
    ft = T(f(-1 + t))
    y0 = K(-a + b - c + 1)
    f1 = ft[1]
    f2 = ft[2]
    f3 = ft[3]
    y1 = f1 / (2 * y0)
    y2 = (f2 - y1**2) / (2 * y0)
    y3 = (f3 - 2 * y1 * y2) / (2 * y0)
    yseries = y0 + y1 * t + y2 * t**2 + y3 * t**3

    for k in range(4):
        row = [K(0)] * 8
        for i in range(5):
            row[i] += T((-1 + t) ** i)[k]
        for j in range(3):
            row[5 + j] += T(((-1 + t) ** j) * yseries)[k]
        rows.append(row)

    # Vanishing at W=(r,0).
    rows.append([K(1), K(r), K(r) ** 2, K(r) ** 3, K(r) ** 4, K(0), K(0), K(0)])

    poly_rows = []
    for row in rows:
        den = R(1)
        for entry in row:
            den = lcm_poly(den, R(entry.denominator()))
        poly_rows.append([R(entry * den) for entry in row])

    dets = []
    for omit in range(9):
        print(f"computing 8x8 minor omitting row {omit}", flush=True)
        sub = Matrix(R, [poly_rows[i] for i in range(9) if i != omit])
        dets.append(R(sub.det()))

    print("computing common gcd", flush=True)
    common = gcd([d for d in dets if d])
    detcores = [R(d // common) for d in dets]

    lines.append("Contact model:")
    lines.append("  h = a*x^3 + b*x^2 + c*x + 1")
    lines.append("  f = h^2 - a^2*x^5*(x+1)")
    lines.append("  div(y-h) = 5*(0,1) + (-1,h(-1)) - 6*infinity")
    lines.append("")
    lines.append("Rank-drop matrix rows:")
    lines.append("  0..3: A+B*y vanishes to order 4 at (0,1)")
    lines.append("  4..7: A+B*y vanishes to order 4 at (-1,h(-1))")
    lines.append("  8:    A(r)=0 at W=(r,0)")
    lines.append("")
    lines.append(f"Common factor in all 8x8 minors: {common.factor()}")
    lines.append("")

    # The minor omitting row 8 is the condition for a function vanishing to
    # order 4 at both P0 and P1.  It is linear after p=a-b+c-1 is introduced.
    H = detcores[8]

    S = PolynomialRing(QQ, ("A", "P", "C", "Rr"))
    A, P, C, Rr = S.gens()
    phi = R.hom([A, A + C - 1 - P, C, Rr], S)
    Hp = phi(H)

    T3 = PolynomialRing(QQ, ("A", "P", "Rr"))
    A3, P3, R3 = T3.gens()
    H_C0 = coeff_by_var(Hp, 2, 0, T3)
    H_C1 = coeff_by_var(Hp, 2, 1, T3)
    H_C2 = coeff_by_var(Hp, 2, 2, T3)
    assert H_C2 == 0
    C_expr = -T3.fraction_field()(H_C0) / T3.fraction_field()(H_C1)

    lines.append("First non-boundary rank condition after p=a-b+c-1:")
    lines.append(f"  H(A,P,C) = {Hp}")
    lines.append("  This is linear in C, with")
    lines.append(f"  C = {C_expr}")
    lines.append("  and b = A + C - 1 - P.")
    lines.append("")

    # Substitute the solved C into f(r)=0 and the remaining rank minors.
    KT3 = T3.fraction_field()
    Aq, Pq, Rq = [KT3(g) for g in (A3, P3, R3)]
    Cq = KT3(C_expr)
    psi = R.hom([Aq, Aq + Cq - 1 - Pq, Cq, Rq], KT3)

    fR = R(f(r))
    fr_sub = T3(psi(fR).numerator())
    minor_subs = [T3(psi(d).numerator()) for d in detcores]

    printed = {
        A3: QQ(-240),
        P3: QQ(-1440) / 11,
        R3: QQ(1) / 15,
    }
    printed_C = C_expr.subs(printed)

    lines.append("Printed Elkies member in this normalization:")
    lines.append("  h_unscaled = -2640*x^3 - 1323*x^2 - 112*x + 11")
    lines.append("  after y -> y/11:")
    lines.append("  A=a=-240, P=p=-1440/11, R=r=1/15")
    lines.append(f"  C from formula = {printed_C}")
    lines.append("  expected C=c=-112/11")
    lines.append(f"  f(R) numerator at printed point = {fr_sub.subs(printed)}")
    lines.append("  remaining minor numerators at printed point:")
    lines.append("    " + ", ".join(str(m.subs(printed)) for m in minor_subs))
    lines.append("")

    # Homogenize by the scale ratio z=A/P.  The printed example has z=11/6.
    U = PolynomialRing(QQ, ("z", "P", "Rr"))
    z, Pu, Ru = U.gens()
    theta = T3.hom([z * Pu, Pu, Ru], U)
    Fz = theta(fr_sub)
    Gz_candidates = []
    for i, m in enumerate(minor_subs):
        if i == 8:
            continue
        gz = theta(m)
        if gz != 0 and gz.subs({z: QQ(11) / 6, Pu: QQ(-1440) / 11, Ru: QQ(1) / 15}) == 0:
            Gz_candidates.append((i, gz))

    lines.append("After setting z=A/P:")
    lines.append(f"  deg_P f(R) numerator = {Fz.degree(Pu)}")
    lines.append("  rank minors that still vanish on the printed point:")
    lines.append("    " + ", ".join(f"omit row {i}, deg_P={g.degree(Pu)}" for i, g in Gz_candidates))
    lines.append("")

    def factors_vanishing_at_printed(poly):
        fac = poly.factor()
        point = {z: QQ(11) / 6, Pu: QQ(-1440) / 11, Ru: QQ(1) / 15}
        return [(base, exp) for base, exp in fac if base.subs(point) == 0]

    F_vanish = factors_vanishing_at_printed(Fz)
    lines.append("Factors of f(R) numerator vanishing at printed point:")
    for base, exp in F_vanish:
        lines.append(f"  exponent {exp}, degrees (z,P,R)=({base.degree(z)},{base.degree(Pu)},{base.degree(Ru)}):")
        lines.append(f"    {base}")
    lines.append("")

    # Choose the first nonzero remaining minor that has a vanishing factor.
    chosen_i = None
    chosen_G = None
    G_vanish = []
    for i, gz in Gz_candidates:
        vanish = factors_vanishing_at_printed(gz)
        if vanish:
            chosen_i = i
            chosen_G = gz
            G_vanish = vanish
            break

    if chosen_i is None:
        lines.append("No remaining rank minor factor vanished at the printed point.")
    else:
        lines.append(f"Using remaining rank minor omitting row {chosen_i}.")
        lines.append("Factors of that minor vanishing at printed point:")
        for base, exp in G_vanish:
            lines.append(f"  exponent {exp}, degrees (z,P,R)=({base.degree(z)},{base.degree(Pu)},{base.degree(Ru)}):")
            lines.append(f"    {base}")
        lines.append("")

        F_component = prod(base**exp for base, exp in F_vanish)
        G_component = prod(base**exp for base, exp in G_vanish)
        F_component = U(F_component)
        G_component = U(G_component)
        lines.append("Component equations before eliminating P:")
        lines.append(f"  F(z,P,R) = {F_component}")
        lines.append(f"  G(z,P,R) = {G_component}")
        lines.append("")

        V = PolynomialRing(QQ, ("z", "Rr"))
        zv, Rv = V.gens()
        F2 = coeff_by_var(F_component, 1, 2, V)
        F1 = coeff_by_var(F_component, 1, 1, V)
        F0 = coeff_by_var(F_component, 1, 0, V)
        G2 = coeff_by_var(G_component, 1, 2, V)
        G1 = coeff_by_var(G_component, 1, 1, V)
        G0 = coeff_by_var(G_component, 1, 0, V)
        P_num = F2 * G0 - G2 * F0
        P_den = G2 * F1 - F2 * G1
        pc = gcd(P_num, P_den)
        if pc != 0:
            P_num //= pc
            P_den //= pc
        printed_zr = {zv: QQ(11) / 6, Rv: QQ(1) / 15}
        lines.append("Common-root formula for the scale P=p:")
        lines.append("  If F=F2*P^2+F1*P+F0 and G=G2*P^2+G1*P+G0, then")
        lines.append("  P = (F2*G0 - G2*F0)/(G2*F1 - F2*G1) on the component.")
        lines.append(f"  numerator factorization: {short_factor_line(P_num.factor())}")
        lines.append(f"  denominator factorization: {short_factor_line(P_den.factor())}")
        lines.append(f"  P at printed point = {P_num.subs(printed_zr) / P_den.subs(printed_zr)}")
        lines.append("")

        res = F_component.resultant(G_component, Pu)
        fac_res = res.factor()
        point_zr = {z: QQ(11) / 6, Ru: QQ(1) / 15}
        lines.append("Resultant Res_P(F,G) factorization summary:")
        for base, exp in fac_res:
            val = base.subs(point_zr)
            flag = "PRINTED_COMPONENT" if val == 0 else "nonzero_at_printed"
            lines.append(
                f"  exponent {exp}, degrees (z,R)=({base.degree(z)},{base.degree(Ru)}), {flag}:"
            )
            lines.append(f"    {base}")
        lines.append("")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
