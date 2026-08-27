#!/usr/bin/env sage --python
"""Geometric simplicity check for the M(2,12)+extra-3 hit."""

from sage.all import GF, HyperellipticCurve, PolynomialRing, QQ, EllipticCurve, prime_range

R = PolynomialRing(QQ, "x")
x = R.gen()

f = (
    5668704*x**5 - 22143375*x**4 + 36098622*x**3
    - 30305259*x**2 + 12990780*x - 2259900
)
C = HyperellipticCurve(f)
J = C.jacobian()

E1 = EllipticCurve("90c1")
E2 = EllipticCurve("510g1")

print("geometric_endomorphism_algebra_is_field(B=100) =",
      J.geometric_endomorphism_algebra_is_field(B=100))
print("geometric_endomorphism_ring_is_ZZ(B=100) =",
      J.geometric_endomorphism_ring_is_ZZ(B=100))
print("E1", E1.cremona_label(), "ainvs", E1.ainvs(), "conductor", E1.conductor())
print("E2", E2.cremona_label(), "ainvs", E2.ainvs(), "conductor", E2.conductor())

bad = {2, 3, 5, 17}
checked = []
for p in prime_range(7, 300):
    if p in bad:
        continue
    char = C.change_ring(GF(p)).frobenius_polynomial()
    target = (x**2 - E1.ap(p)*x + p) * (x**2 - E2.ap(p)*x + p)
    if char != target:
        raise RuntimeError(f"Frobenius mismatch at p={p}: {char} != {target}")
    checked.append(p)

print("matched product Frobenius factors through p =", checked[-1])
print("number of primes checked =", len(checked))
