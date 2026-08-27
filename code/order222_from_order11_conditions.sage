from sage.all import *
import sys

# Extra 2-torsion test on the order-22 subfamily obtained by forcing
# a rational branch point in Flynn's order-11 infinity-torsion family.

R = PolynomialRing(QQ, ("s", "x"))
s, x = R.gens()

Px = PolynomialRing(QQ, "x")
X = Px.gen()


def q_minus(S=s, Xvar=x):
    return (-(S - 1) ** 2 * Xvar ** 2
            + (S - 1) * (S ** 3 - S ** 2 + S + 1) * Xvar
            - (S ** 2 - S + 1) ** 2)


def q_plus(S=s, Xvar=x):
    return (-(S + 1) ** 2 * Xvar ** 2
            + (S + 1) * (S ** 3 + S ** 2 + S - 1) * Xvar
            - (S ** 2 + S + 1) ** 2)


def c_minus(S=s, Xvar=x):
    return (
        S ** 6 * Xvar - S ** 6 - 2 * S ** 5 * Xvar
        + 2 * S ** 4 * Xvar ** 2 + 2 * S ** 5 + 3 * S ** 4 * Xvar
        - 4 * S ** 3 * Xvar ** 2 + S ** 2 * Xvar ** 3
        - 3 * S ** 4 - 2 * S ** 3 * Xvar + 5 * S ** 2 * Xvar ** 2
        - 2 * S * Xvar ** 3 + 2 * S ** 3 + 2 * S ** 2 * Xvar
        - 4 * S * Xvar ** 2 + Xvar ** 3 - S ** 2 - 2 * S * Xvar
        + Xvar ** 2 + Xvar
    )


def c_plus(S=s, Xvar=x):
    return c_minus(-S, Xvar)


def qdisc(eps, S=s):
    if eps == -1:
        return (S - 1) ** 5 * (S ** 3 + S ** 2 - S + 3)
    if eps == 1:
        return (S + 1) ** 5 * (S ** 3 - S ** 2 - S - 3)
    raise ValueError("eps must be +/-1")


def qdisc_quartic(eps, S=s):
    # Away from the boundary S=-eps, the fifth power differs from this
    # quartic condition by a square factor.
    return S ** 4 - 2 * S ** 2 - 4 * eps * S - 3


def rat_height_values(B):
    for den in range(1, B + 1):
        for num in range(-B, B + 1):
            if gcd(num, den) == 1 and max(abs(num), den) <= B:
                yield QQ(num) / QQ(den)


def is_square_qq(a):
    a = QQ(a)
    return a >= 0 and ZZ(a.numerator()).is_square() and ZZ(a.denominator()).is_square()


def cubic_poly(eps, S):
    S = QQ(S)
    if eps == -1:
        return c_minus(S, X)
    if eps == 1:
        return c_plus(S, X)
    raise ValueError("eps must be +/-1")


def rational_roots(poly):
    roots = []
    for factor, mult in poly.factor():
        if factor.degree() == 1:
            roots.append(-factor[0] / factor[1])
    return roots


def residual_search(B):
    hits = []
    for S in rat_height_values(B):
        if S == 0 or S ** 2 == 1:
            continue
        for eps in [-1, 1]:
            D = qdisc(eps, S)
            if is_square_qq(D):
                hits.append((eps, S, "quadratic_split", D))
            roots = rational_roots(cubic_poly(eps, S))
            if roots:
                hits.append((eps, S, "cubic_root", roots))
    return hits


def report_geometry():
    A2 = AffineSpace(QQ, 2, names=("s", "x"))
    ss, xx = A2.gens()
    for eps, cfun in [(-1, c_minus), (1, c_plus)]:
        Caff = Curve(cfun(ss, xx))
        Cproj = Caff.projective_closure()
        print("cubic residual eps", eps)
        print("  projective genus", Cproj.genus())
        print("  projective closure singular", Cproj.is_singular())


print("Residual quadratic factors:")
print("q_minus =", q_minus())
print("disc_minus =", qdisc(-1))
print("square condition minus: Y^2 =", qdisc_quartic(-1))
print("q_plus =", q_plus())
print("disc_plus =", qdisc(1))
print("square condition plus: Y^2 =", qdisc_quartic(1))
print()
report_geometry()
print()

B = Integer(sys.argv[1]) if len(sys.argv) > 1 else Integer(500)
hits = residual_search(B)
print("B", B, "hits", hits, "count", len(hits))
