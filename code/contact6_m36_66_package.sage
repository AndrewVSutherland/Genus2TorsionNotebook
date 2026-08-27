#!/usr/bin/env sage
"""Exact certificate for the contact-6 [6,6] specialization.

Run with:
    sage code/contact6_m36_66_package.sage
"""

from sage.all import GF, GCD, HyperellipticCurve, PolynomialRing, QQ, ZZ, pari


R = PolynomialRing(QQ, "x")
x = R.gen()

a = QQ(133) / 39
b = -QQ(7) / 13
L = QQ(29) / 16
U = -QQ(9) / 4
v0 = QQ(5) / 2

h6 = 1 + a*x + b*x**2 + x**3
f = h6**2 - (x - 1)**6
den = 1521
F_integral = R(den**2 * f)

# PARI minimisation divides the square factor 78^2 from this integral model.
minimal_data = pari(F_integral).hyperellminimalmodel()
F = R(F_integral / 78**2)
assert minimal_data[0] == pari(F)
assert minimal_data[1] == 0
assert pari(F).hyperellred()[0] == pari(F)
assert F == 1872*x**5 - 3000*x**4 + 6969*x**3 - 1691*x**2 + 4875*x

scale = QQ(39) / 2
assert F == scale**2 * f

# Recover the cubic-contact class and transfer both contacts to the minimal model.
M = L**2
c = [f[i] for i in range(6)]
B = c[5]*M + 3*U
Delta = 4*c[4]*M + 12*(U**2 + v0**2) - B**2
q = x**2 + U*x + v0**2
h3 = (x**3 + (B/2)*x**2 + (Delta/8)*x + v0**3) / L
H6 = scale*h6
H3 = scale*h3

assert H6**2 - F == scale**2*(x - 1)**6
assert H3**2 - F == (QQ(97344)/841)*q**3
assert H3 % q == 25*x - 300

C = HyperellipticCurve(F)
JQ = C.jacobian()(QQ)
zero = JQ(0)

D = JQ([x - 1, R(95)])
E = JQ([q, 25*x - 300])
T = JQ([x, R(0)])
W = E + T


def point_order(P, bound=12):
    for n in range(1, bound + 1):
        if n*P == zero:
            return n
    raise RuntimeError("point order exceeds the supplied bound")


assert point_order(D) == 6
assert point_order(E) == 3
assert point_order(T) == 2
assert point_order(W) == 6

subgroup = []
for i in range(6):
    for j in range(6):
        point = i*D + j*W
        if not any(point == old for old in subgroup):
            subgroup.append(point)
assert len(subgroup) == 36

# At good reduction, rational torsion injects away from the residue
# characteristic.  The two counts bound the full torsion order by their gcd.
reduction_data = []
for p in (7, 11):
    assert ZZ(F.discriminant()) % p != 0
    fp = C.change_ring(GF(p)).frobenius_polynomial()
    reduction_data.append((p, fp, ZZ(fp(1))))
assert GCD([row[2] for row in reduction_data]) == 36

# Reproduce the 12th-power Frobenius transform in Lombardo's Algorithm 4.10.
R2 = PolynomialRing(ZZ, 2, "xv")
xx, vv = R2.gens()
Tv = PolynomialRing(QQ, "v")
fp37 = C.change_ring(GF(37)).frobenius_polynomial()
fp37_12 = Tv(R2(fp37).resultant(vv - xx**12))
assert fp37_12.is_irreducible()

J = C.jacobian()
geom_field = J.geometric_endomorphism_algebra_is_field(B=100)
geom_ZZ = J.geometric_endomorphism_ring_is_ZZ(B=100)
assert geom_field and geom_ZZ

print("contact parameters:", "a =", a, "b =", b, "L =", L, "U =", U, "v =", v0)
print("earlier square-scaled integral model:", F_integral)
print("PARI minimal and reduced model:", F)
print("minimal model factorization:", F.factor())
print("minimal model discriminant factorization:", ZZ(F.discriminant()).factor())
print("marked D:", D, "order", point_order(D))
print("marked E:", E, "order", point_order(E))
print("rational T:", T, "order", point_order(T))
print("second generator W = E + T:", W, "order", point_order(W))
print("number of distinct i*D+j*W, 0 <= i,j < 6:", len(subgroup))
for p, fp, count in reduction_data:
    print("p =", p, "Frobenius =", fp, "#J(F_p) =", count)
print("gcd of reduction orders:", GCD([row[2] for row in reduction_data]))
print("therefore J(Q)_torsion invariants: [6, 6]")
print("Igusa-Clebsch invariants:", C.igusa_clebsch_invariants())
print("absolute Igusa invariants (Kohel):", C.absolute_igusa_invariants_kohel())
print("p = 37 Frobenius:", fp37)
print("p = 37 12th-power transform:", fp37_12)
print("p = 37 12th-power transform irreducible:", fp37_12.is_irreducible())
print("geometric_endomorphism_algebra_is_field(B=100):", geom_field)
print("geometric_endomorphism_ring_is_ZZ(B=100):", geom_ZZ)
