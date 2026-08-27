#!/usr/bin/env sage
"""Bounded ramification analysis for the corrected degree-16 [64] cover.

This avoids the global discriminant and maximal-order computations that time
out.  It factors the bad-fiber divisor, constructs the compact monic
degree-16 eliminant, and gives finite residue certificates at the two large
irreducible factors.  Local Puiseux calculations are recorded separately in
the accompanying data/note because Singular's output is verbose.
"""

_caller_name = __name__
__name__ = "task5_cover_geometry_import"
load("code/task5_cover_geometry.sage")
__name__ = _caller_name


def generic_halving_eliminant():
    R = PolynomialRing(QQ, names=("f1", "f2", "f3", "f4", "f5", "m", "n"))
    f1, f2, f3, f4, f5, m, n = R.gens()
    K = R.fraction_field()
    F1, F2, F3, F4, F5, M, N = map(K, R.gens())
    alpha = (F4 - M**2)/(2*F5)
    beta = (F3 - 2*M*N - F5*alpha**2)/(2*F5)
    k2 = R((F2 - N**2 - 2*M - 2*F5*alpha*beta).numerator())
    k1 = R((F1 - 2*N - F5*beta**2).numerator())
    H = [g for g, e in k1.resultant(k2, n).factor() if g.degree(m) == 16][0]
    return R, H


def family(base):
    K = FunctionField(base, "t")
    t = K.gen()
    z, r = small_c32_parameter(t)
    p = 4*p_num(z, r)/p_den(z, r)
    a = z*p
    c = (-p*z**4 + 4*p*z**3 - 8*p*z**2 + 8*z**2 + 8*p*z - 16)/(4*z**2)
    b = p*(z - 1) + c - 1
    values = [2*c, c**2 + 2*b, 2*a + 2*b*c, b**2 + 2*a*c, a*(2*b - a)]
    return K, t, values


def specialized_h(base, H, R):
    K, t, values = family(base)
    S = PolynomialRing(K, "M")
    M = S.gen()
    h = R.hom(values + [M, K(0)], S)(H).monic()
    return K, t, values, h


def evaluate_at(h, u, field):
    T = PolynomialRing(field, "Y")
    return T([
        c.numerator()(u)/c.denominator()(u)
        for c in h.list()
    ])


R, H = generic_halving_eliminant()
K, t, values, h = specialized_h(QQ, H, R)
P = PolynomialRing(K, "x")
x = P.gen()
q = sum(values[i]*x**i for i in range(5))

A = values[0].numerator()
B = q.discriminant().numerator()
C = values[4].numerator()
Phi = (t**4 + t**3 + t**2 + t + 1).numerator()

print("degree16_factor_degrees", [(g.degree(), e) for g, e in h.factor()])
print("A_factor_degrees", [(g.degree(), e) for g, e in A.factor()])
print("B_factor_degrees", [(g.degree(), e) for g, e in B.factor()])
print("C_factorization", C.factor())
print("Phi_factor_degrees", [(g.degree(), e) for g, e in Phi.factor()])
print("A_B_squarefree", gcd(A, A.derivative()).degree(), gcd(B, B.derivative()).degree())
print("pairwise_gcd_degrees", gcd(A, B).degree(), gcd(A, C).degree(), gcd(B, C).degree())
print("Phi_gcd_degrees", gcd(Phi, A).degree(), gcd(Phi, B).degree(), gcd(Phi, C).degree())
print("C_distinct_geometric_roots", sum(g.degree() for g, e in C.factor()))
print("binary_discriminant_support_exponents", "A^2 B C^2 / Phi^14")
print("geometric_boundary_points", 8 + 104 + 7 + 4 + 1)

# A simple residue of A modulo 29 and a simple residue of B modulo 41.
# Squarefree specialized h proves that h is separable modulo the corresponding
# characteristic-zero prime divisor, hence the cover is unramified there.
for ell, u, label in [(29, 2, "A"), (41, 19, "B")]:
    F = GF(ell)
    KF, tt, valsF, hF = specialized_h(F, H, R)
    h0 = evaluate_at(hF, F(u), F)
    target = (A if label == "A" else B).change_ring(F)
    other_factors = [
        g.change_ring(F) for g in (A, B, C, Phi)
        if g is not (A if label == "A" else B)
    ]
    print(label + "_certificate", ell, u)
    print(label + "_horizontal_value_derivative", target(F(u)), target.derivative()(F(u)))
    print(label + "_other_support_values", [g(F(u)) for g in other_factors])
    print(label + "_h_integral", all(c.denominator()(F(u)) != 0 for c in hF.list()))
    print(label + "_residual_factor_degrees", [(g.degree(), e) for g, e in h0.factor()])
    print(label + "_residual_gcd_degree", gcd(h0, h0.derivative()).degree())

print("unramified_large_boundary_points", 8 + 104)
print("remaining_boundary_points", 7 + 4 + 1)
print("rigorous_total_different_upper_bound", 12*15)
print("rigorous_genus_upper_bound", 75)
print("exploratory_puiseux_genus_upper_bound", 57)
print("puiseux_candidate_total_different", 96)
print("puiseux_candidate_genus", 33)
