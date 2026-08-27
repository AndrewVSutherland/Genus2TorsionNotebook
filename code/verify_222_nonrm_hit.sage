#!/usr/bin/env sage
"""[2,22] hunt: certified End(Jac_Qbar) = ZZ check for a sieve hit on
f = (x^2-1)(x^2+A x+B)(x^2+C x+D).  This is the decisive NON-RM
certificate (geometric-simplicity Frobenius certs cannot distinguish
End = ZZ from real multiplication).

Run: /opt/sage/envs/sage/bin/sage code/verify_222_nonrm_hit.sage A B C D
"""
import sys
from sage.all import HyperellipticCurve, PolynomialRing, QQ

A, B, C, D = (int(t) for t in sys.argv[1:5])
R = PolynomialRing(QQ, "x")
x = R.gen()
f = (x**2 - 1) * (x**2 + A*x + B) * (x**2 + C*x + D)
print("tuple:", (A, B, C, D))
print("f =", f)
Ccrv = HyperellipticCurve(f)
J = Ccrv.jacobian()

alg_is_field = J.geometric_endomorphism_algebra_is_field(B=100)
ring_is_ZZ = J.geometric_endomorphism_ring_is_ZZ(B=100)
print("geometric_endomorphism_algebra_is_field(B=100) =", alg_is_field)
print("geometric_endomorphism_ring_is_ZZ(B=100) =", ring_is_ZZ)

if not ring_is_ZZ:
    raise RuntimeError("FAIL: End(Jac_Qbar) not certified ZZ -- possibly RM")
print("NONRM_CERTIFIED: End(Jac_Qbar) = ZZ (geometrically simple, no RM)")
