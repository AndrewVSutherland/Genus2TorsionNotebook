# Geometric-simplicity check for the t=1 member of Leprévost's
# one-parameter order-21 family.

from sage.all import *

R.<x> = PolynomialRing(QQ)
f = -216*x^5 + 657*x^4 - 696*x^3 + 466*x^2 - 224*x + 49
assert f.discriminant() != 0

C = HyperellipticCurve(f)
J = C.jacobian()
print("geometric_endomorphism_algebra_is_field(B=100) =",
      J.geometric_endomorphism_algebra_is_field(B=100))
print("geometric_endomorphism_ring_is_ZZ(B=100) =",
      J.geometric_endomorphism_ring_is_ZZ(B=100))
