#!/usr/bin/env sage -python
"""Algebraic extra 2- and 3-torsion conditions for the order-30 family."""

from __future__ import annotations

from math import gcd

from sage.all import QQ, PolynomialRing, FractionField, factor


def is_square_q(q) -> bool:
    q = QQ(q)
    return q >= 0 and q.numerator().is_square() and q.denominator().is_square()


def rational_params(bound: int):
    seen = set()
    for den in range(1, bound + 1):
        for num in range(-bound, bound + 1):
            if gcd(num, den) != 1:
                continue
            val = QQ(num) / QQ(den)
            if val not in seen:
                seen.add(val)
                yield val


def family_data_numeric(Rv, branch):
    P = PolynomialRing(QQ, "x")
    x = P.gen()
    if Rv**2 == 5:
        return None
    t = (5 * Rv**2 - 20 * Rv + 19) / (Rv**2 - 5)
    Y = -2 * (5 * Rv**2 - 22 * Rv + 25) / (Rv**2 - 5)
    u = t**3
    if u == 0:
        return None
    s = (
        t**5
        + t**4
        + QQ(5) / 2 * t**3
        + QQ(1) / 2 * t
        + branch * t * (t - QQ(1) / 2) * (t + 1) * Y
    )
    C = (u**2 + 1) / (2 * u)
    c = (u**2 - 1) / (2 * u)
    denq = u**6 + 6 * u**4 * s - 2 * u**4 + 15 * u**3 * s - u * s**3 + u**2
    if c == 0 or denq == 0:
        return None
    numq = (
        15 * u**5
        + 90 * u**4
        + 20 * u**3 * s
        - 6 * u**2 * s**2
        + 231 * u**3
        + 2 * u**2 * s
        - 15 * u * s**2
        + 90 * u**2
        - 20 * u * s
        + 15 * u
        - 2 * s
    )
    q = numq / denq
    A = (s + q) / 2
    B = (15 - s * q) / 2
    D2 = (B - 3) ** 2 - 4 * (A + 3) * (C + 1)
    cubic = 2 * x**3 + (A - 3) * x**2 + (B + 3) * x + (C - 1)
    return t, D2, cubic


def main() -> None:
    Rring = PolynomialRing(QQ, ("R", "rho"))
    Rv, rho = Rring.gens()
    K = FractionField(Rring)

    lines: list[str] = []
    lines.append("# Extra torsion conditions for contact5/contact6 order-30 family")
    lines.append("")
    lines.append("For every specialization,")
    lines.append("  f = Q2*C3")
    lines.append("with")
    lines.append("  Q2 = (A+3)*x^2 + (B-3)*x + (C+1),")
    lines.append("  C3 = 2*x^3 + (A-3)*x^2 + (B+3)*x + (C-1).")
    lines.append("The built-in rational 2-torsion comes from this 2+3 factorization.")
    lines.append("An additional independent 2-torsion class appears iff")
    lines.append("  (i) Q2 splits: Delta2 is a square, or")
    lines.append("  (ii) C3 has a rational root rho.")
    lines.append("")

    for eps in [-1, 1]:
        r = Rv
        t = (5 * r**2 - 20 * r + 19) / (r**2 - 5)
        Y = -2 * (5 * r**2 - 22 * r + 25) / (r**2 - 5)
        u = t**3
        s = (
            t**5
            + t**4
            + QQ(5) / 2 * t**3
            + QQ(1) / 2 * t
            + eps * t * (t - QQ(1) / 2) * (t + 1) * Y
        )
        C = (u**2 + 1) / (2 * u)
        c = (u**2 - 1) / (2 * u)
        denq = u**6 + 6 * u**4 * s - 2 * u**4 + 15 * u**3 * s - u * s**3 + u**2
        numq = (
            15 * u**5
            + 90 * u**4
            + 20 * u**3 * s
            - 6 * u**2 * s**2
            + 231 * u**3
            + 2 * u**2 * s
            - 15 * u * s**2
            + 90 * u**2
            - 20 * u * s
            + 15 * u
            - 2 * s
        )
        q = numq / denq
        A = (s + q) / 2
        B = (15 - s * q) / 2

        delta2 = K((B - 3) ** 2 - 4 * (A + 3) * (C + 1))
        delta_num = Rring(delta2.numerator())
        delta_den = Rring(delta2.denominator())
        c3_root = K(2 * rho**3 + (A - 3) * rho**2 + (B + 3) * rho + (C - 1))
        c3_num = Rring(c3_root.numerator())

        lines.append(f"BRANCH eps={eps}")
        lines.append("Delta2 numerator factorization:")
        lines.append(f"  {factor(delta_num)}")
        lines.append("Delta2 denominator factorization:")
        lines.append(f"  {factor(delta_den)}")
        lines.append("Reduced squareclass condition for Q2 splitting:")
        if eps == -1:
            lines.append(
                "  W^2 = -1095*(R^2 - 4*R + 13/3)"
                "*(R^6 - 4044/365*R^5 + 18249/365*R^4"
                " - 42664/365*R^3 + 54039/365*R^2"
                " - 34764/365*R + 1751/73)"
            )
        else:
            lines.append(
                "  W^2 = -2135*(R^2 - 32/7*R + 37/7)"
                "*(R^6 - 3504/305*R^5 + 16761/305*R^4"
                " - 42784/305*R^3 + 61611/305*R^2"
                " - 47664/305*R + 3119/61)"
            )
        lines.append(
            "C3 rational-root condition: numerator of "
            "2*rho^3+(A-3)*rho^2+(B+3)*rho+(C-1) is:"
        )
        lines.append(f"  degree_R={c3_num.degree(Rv)} degree_rho={c3_num.degree(rho)}")
        lines.append(f"  factorization={factor(c3_num)}")
        lines.append("")

    lines.append("Additional independent 3-torsion condition")
    lines.append("Let f=sum c_i*x^i and seek q3=x^2+U*x+V,")
    lines.append("ell3=M*x^3+N*x^2+P*x+S with")
    lines.append("  ell3^2 - f = M^2*q3^3.")
    lines.append("The six coefficient equations are:")
    lines.append("  2*M*N - c5 - 3*M^2*U = 0")
    lines.append("  N^2 + 2*M*P - c4 - M^2*(3*U^2+3*V) = 0")
    lines.append("  2*M*S + 2*N*P - c3 - M^2*(U^3+6*U*V) = 0")
    lines.append("  P^2 + 2*N*S - c2 - M^2*(3*U^2*V+3*V^2) = 0")
    lines.append("  2*P*S - c1 - 3*M^2*U*V^2 = 0")
    lines.append("  S^2 - c0 - M^2*V^3 = 0")
    lines.append("For this order-30 family, substitute the c_i from")
    lines.append("  f=(x^3+A*x^2+B*x+C)^2-(x-1)^6.")
    lines.append("The known 3-torsion is the non-reduced solution")
    lines.append("  q3=(x-1)^2, i.e. U=-2, V=1, ell3=+/-h6.")
    lines.append("A genuinely additional 3-torsion class must avoid that")
    lines.append("branch and pass an independence check in J[3].")
    lines.append("")

    height = 120
    hits = []
    for Rtest in rational_params(height):
        for eps in [-1, 1]:
            dat = family_data_numeric(Rtest, eps)
            if dat is None:
                continue
            tval, delta2, cubic = dat
            if is_square_q(delta2):
                hits.append(("Q2split", Rtest, eps, tval, delta2, None))
            roots = cubic.roots(QQ)
            if roots:
                hits.append(("C3root", Rtest, eps, tval, None, roots))
    lines.append(f"Small exact extra-2 scan: height {height}, hits {len(hits)}.")
    for hit in hits[:50]:
        lines.append(f"  {hit}")

    with open("data/contact5_contact6_order30_extra_conditions.txt", "w", encoding="ascii") as fh:
        fh.write("\n".join(lines) + "\n")
    print("Wrote data/contact5_contact6_order30_extra_conditions.txt")


if __name__ == "__main__":
    main()
