# Sage/Lombardo geometric simplicity checks for published order-11 families.

from sage.all import *

R.<x> = PolynomialRing(QQ)


def flynn11(t):
    t = QQ(t)
    return x**6 + 2*x**5 + (2*t + 3)*x**4 + 2*x**3 + (t**2 + 1)*x**2 + 2*t*(1 - t)*x + t**2


def daowsud_schmidt11(u):
    u = QQ(u)
    return (x**6 - 4*x**5 + 8*(1 + u)*x**4 - (10 + 32*u)*x**3
            + 8*(1 + 6*u + 2*u**2)*x**2
            - 4*(1 + 6*u + 16*u**2)*x + 64*u**2 + 1)


def check(label, f):
    print(label)
    print("disc_zero", f.discriminant() == 0)
    C = HyperellipticCurve(f)
    J = C.jacobian()
    print("geometric_endomorphism_algebra_is_field(B=100) =", J.geometric_endomorphism_algebra_is_field(B=100))
    print("geometric_endomorphism_ring_is_ZZ(B=100) =", J.geometric_endomorphism_ring_is_ZZ(B=100))


for t in [1, 2, QQ(1)/2, 3]:
    check(f"Flynn t={t}", flynn11(t))

for u in [1, 2, -1, QQ(1)/2, 3]:
    check(f"Daowsud-Schmidt u={u}", daowsud_schmidt11(u))
