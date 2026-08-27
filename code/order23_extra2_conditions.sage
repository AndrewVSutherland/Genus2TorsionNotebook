# Derive the rational-2-torsion condition on the Kuru-Sadek [23] family.
# A finite rational branch point r is equivalent to a rational root of f_t(x).
# Since the family is invariant under t -> -t, the condition descends to q=t^2.

from sage.all import *

T.<u> = PolynomialRing(QQ)
K = FractionField(T)
R.<x> = PolynomialRing(K)
S.<q, r> = PolynomialRing(QQ, 2)


def polynomial_to_q_ring(poly):
    out = S(0)
    for j, cj in enumerate(T(poly).list()):
        if cj:
            assert j % 2 == 0
            out += QQ(cj)*q**(j//2)
    return out


def primitive_part(P):
    lcmd = LCM([QQ(c).denominator() for c in P.coefficients()])
    P = S(lcmd)*P
    cont = GCD([ZZ(c) for c in P.coefficients()])
    return P/S(cont)


def kuru_sadek23_function_field():
    beta = (u**2 + 1)**2/(4*u**2)
    sbeta = (u**2 + 1)/(2*u)
    s = (u**2 - 1)/(2*u)
    alpha = K(beta) - K(s**5)/K(beta*sbeta)
    lam = (alpha - 1)**4/((alpha - K(beta))**2*alpha)
    expr = (x**3*(x - alpha)**2
            - (x - 1)*((x - 1)**4 - lam*(x - K(beta))**2*x))/(2*(x - alpha)*(x - K(beta)))
    A, rem = R(expr.numerator()).quo_rem(R(expr.denominator()))
    assert rem == 0
    return R(A**2 - lam*x**4*(x - 1))


def rational_root_condition():
    f = kuru_sadek23_function_field()
    lcmden = T(1)
    for c in f.coefficients(sparse=False):
        lcmden = lcm(lcmden, T(c.denominator()))

    G = S(0)
    for i, c in enumerate(f.coefficients(sparse=False)):
        G += polynomial_to_q_ring(lcmden*c)*r**i
    return primitive_part(G), lcmden


def boundary_q(q0, p):
    q0 = GF(p)(q0)
    if q0 == 0 or q0 == 1 or q0 == -1:
        return True
    if p != 5 and 5*q0**4 + 10*q0**2 + 1 == 0:
        return True
    return False


G, lcmden_t = rational_root_condition()
print('Kuru-Sadek [23] extra rational 2-torsion condition')
print('Let q=t^2 and r be the candidate finite branch point.')
print('A rational branch point is equivalent to G(q,r)=0 with q a rational square.')
print('deg_q', G.degree(q), 'deg_r', G.degree(r), 'total_degree', G.total_degree(), 'terms', len(G.dict()))
print('factorization over QQ:')
print(G.factor())
print('denominator in t:')
print(lcmden_t.factor())
print('coefficient factors by r-power:')
for i in range(G.degree(r) + 1):
    c = G.coefficient({r: i})
    if c:
        print('r^%s:' % i, c.factor())

A2 = AffineSpace(QQ, 2, names=('q', 'r'))
C = Curve(G, A2)
print('affine singular points:')
print(C.singular_points())
print('geometric genus:')
print(C.geometric_genus())

primes = [3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97]
print('finite-prime q-residue filter:')
for p in primes:
    F = GF(p)
    Gp = G.change_ring(F)
    squares = {F(a)**2 for a in F}
    good = []
    pairs = 0
    for q0 in F:
        if q0 not in squares:
            continue
        if boundary_q(q0, p):
            continue
        rs = []
        for r0 in F:
            if Gp(q0, r0) == 0:
                rs.append(ZZ(r0))
        if rs:
            good.append(ZZ(q0))
            pairs += len(rs)
    print('p=%s good_q=%s pairs=%s' % (p, good, pairs))
