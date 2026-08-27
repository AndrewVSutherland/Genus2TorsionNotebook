# Alternative square-root/sign model for rational 2-torsion on the Kuru-Sadek [23] family.
# The root equation f_t(r)=0 implies alpha*(r-1) is a square because lambda*alpha is a square.

from sage.all import *

Qqz = PolynomialRing(QQ, ['q', 'z']); q, z = Qqz.gens()
T = PolynomialRing(QQ, 't'); t = T.gen()
K = FractionField(T)
R = PolynomialRing(K, 'x'); x = R.gen()
Z = PolynomialRing(K, 'zz'); zz = Z.gen()

beta = (t**2 + 1)**2/(4*t**2)
sbeta = (t**2 + 1)/(2*t)
s = (t**2 - 1)/(2*t)
alpha = K(beta) - K(s**5)/(K(beta)*K(sbeta))
lam = (alpha - 1)**4/((alpha - K(beta))**2*alpha)

expr = (x**3*(x-alpha)**2 - (x-1)*((x-1)**4 - lam*(x-K(beta))**2*x))/(2*(x-alpha)*(x-K(beta)))
A, rem = R(expr.numerator()).quo_rem(R(expr.denominator()))
assert rem == 0

def eval_A(val):
    return sum(K(A[i])*val**i for i in range(A.degree()+1))

print('lambda*alpha factor:')
print(K(lam*alpha).factor())
print('alpha factor:')
print(K(alpha).factor())
print('alpha-1 factor:')
print(K(alpha-1).factor())
print('alpha-beta factor:')
print(K(alpha-K(beta)).factor())

# If z^2 = alpha*(r-1), then r = 1 + z^2/alpha.
# The sign equation is A(r)*(alpha-beta)*alpha = (alpha-1)^2*r^2*z.
r0 = 1 + zz**2/alpha
E = Z(eval_A(r0)*(alpha-K(beta))*alpha - (alpha-1)**2*r0**2*zz)
num = E.numerator()
lcmden = T(1)
for c in num.coefficients(sparse=False):
    lcmden = lcm(lcmden, T(c.denominator()))

P = Qqz(0)
for iz, c in enumerate(num.coefficients(sparse=False)):
    poly = T(lcmden*c)
    for it, ct in enumerate(poly.list()):
        if ct:
            assert it % 2 == 0
            P += QQ(ct)*q**(it//2)*z**iz
lcmd = LCM([QQ(c).denominator() for c in P.coefficients()])
P = Qqz(lcmd)*P
cont = GCD([ZZ(c) for c in P.coefficients()])
P = P/Qqz(cont)

print('raw sign equation factor degrees:')
print([(f.degree(q), f.degree(z), e) for f, e in P.factor()])
main = [f for f, e in P.factor() if f.degree(z) > 0][0]
print('main sign model degrees:')
print('deg_q', main.degree(q), 'deg_z', main.degree(z), 'total_degree', main.total_degree(), 'terms', len(main.dict()))

A2 = AffineSpace(QQ, 2, names=('qq', 'zz'))
CR = A2.coordinate_ring(); qq, zz2 = CR.gens()
MM = CR(0)
for (iq, iz), c in main.dict().items():
    MM += c*qq**iq*zz2**iz
C = Curve(MM, A2)
print('affine singular points:')
print(C.singular_points())
print('geometric genus:')
print(C.geometric_genus())

print('simple symmetry tests:')
print('z -> -z proportional:', main(q, -z) == main or main(q, -z) == -main)

def transform_q(poly, num, den):
    deg = poly.degree(q)
    out = Qqz(0)
    for (iq, iz), c in poly.dict().items():
        out += c * num**iq * den**(deg-iq) * z**iz
    return out / gcd(out.coefficients())

for name, num, den in [
    ('1/q', Qqz(1), q),
    ('-q', -q, Qqz(1)),
    ('1-q', 1-q, Qqz(1)),
    ('q/(q-1)', q, q-1),
    ('1/(1-q)', Qqz(1), 1-q),
    ('(q-1)/q', q-1, q),
]:
    Tpoly = transform_q(main, num, den)
    print(name, 'proportional', Tpoly == main or Tpoly == -main, 'gcd_total_degree', gcd(main, Tpoly).total_degree())
