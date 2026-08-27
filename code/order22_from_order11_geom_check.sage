# Sage/Lombardo geometric simplicity checks for the order-22 families
# obtained by forcing a rational branch point in the order-11 families.

from sage.all import *

R.<x> = PolynomialRing(QQ)


def flynn11(t):
    t = QQ(t)
    return x**6 + 2*x**5 + (2*t + 3)*x**4 + 2*x**3 + (t**2 + 1)*x**2 + 2*t*(1 - t)*x + t**2


def ds11(u):
    u = QQ(u)
    return (x**6 - 4*x**5 + 8*(1 + u)*x**4 - (10 + 32*u)*x**3
            + 8*(1 + 6*u + 2*u**2)*x**2
            - 4*(1 + 6*u + 16*u**2)*x + 64*u**2 + 1)


def flynn22_parameter(s, eps):
    s = QQ(s)
    return (-s**2*(s**2 + 1)*(s**4 - s**2 + 1) + 2*eps*s**5)/(s**2 - 1)**2


def ds22_parameter(s, eps):
    s = QQ(s)
    return (-s**2*(s**2 + 1)*(s**4 - s**2 + 1) + 2*eps*s**5)/(4*(s**2 - 1)**2)


def primitive_integral(f):
    coeffs = f.coefficients(sparse=False)
    lcm = LCM([c.denominator() for c in coeffs])
    F = R(lcm*f)
    g = GCD([ZZ(c) for c in F.coefficients(sparse=False)])
    return R(F/g)


def check(label, f):
    F = primitive_integral(f)
    print(label)
    print('primitive', F)
    print('disc_zero', F.discriminant() == 0)
    C = HyperellipticCurve(F)
    J = C.jacobian()
    print('geometric_endomorphism_algebra_is_field(B=100) =', J.geometric_endomorphism_algebra_is_field(B=100))
    print('geometric_endomorphism_ring_is_ZZ(B=100) =', J.geometric_endomorphism_ring_is_ZZ(B=100))


for eps in [-1, 1]:
    s = QQ(2)
    check(f'Flynn22 s={s} eps={eps}', flynn11(flynn22_parameter(s, eps)))

for eps in [-1, 1]:
    s = QQ(2)
    check(f'DaowsudSchmidt22 s={s} eps={eps}', ds11(ds22_parameter(s, eps)))
