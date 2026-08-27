# Sage/Lombardo geometric simplicity check for the Kuru-Sadek order-23 family.

from sage.all import *

R.<x> = PolynomialRing(QQ)


def kuru_sadek23(t):
    t = QQ(t)
    beta = (t**2 + 1)**2/(4*t**2)
    sbeta = (t**2 + 1)/(2*t)
    s = (t**2 - 1)/(2*t)
    alpha = beta - s**5/(beta*sbeta)
    lam = (alpha - 1)**4/((alpha - beta)**2*alpha)
    expr = (x**3*(x - alpha)**2 - (x - 1)*((x - 1)**4 - lam*(x - beta)**2*x))/(2*(x - alpha)*(x - beta))
    num = R(expr.numerator())
    den = R(expr.denominator())
    A, rem = num.quo_rem(den)
    assert rem == 0
    return R(A**2 - lam*x**4*(x - 1))


def primitive_integral(f):
    coeffs = f.coefficients(sparse=False)
    lcm = LCM([c.denominator() for c in coeffs])
    F = R(lcm*f)
    g = GCD([ZZ(c) for c in F.coefficients(sparse=False)])
    return R(F/g)


f = kuru_sadek23(2)
print("primitive integral model:")
print(primitive_integral(f))
print("disc_zero", f.discriminant() == 0)
C = HyperellipticCurve(f)
J = C.jacobian()
print("geometric_endomorphism_algebra_is_field(B=100) =", J.geometric_endomorphism_algebra_is_field(B=100))
print("geometric_endomorphism_ring_is_ZZ(B=100) =", J.geometric_endomorphism_ring_is_ZZ(B=100))
