#!/usr/bin/env sage
"""Sage/Lombardo endomorphism cross-check for (2,2,2,12) curve #3.

Independent of the Magma Frobenius root-power certificates: certifies
End(Jac(C)_Qbar) via Sage's certified genus-2 endomorphism routines
(the same lane as code/m212_extra3_geom_simple_check.py).

Curve: recorded reduced minimal model from data/22212_abc_curve3.txt,
passed as the integral simplified model y^2 = 4f + h^2.

Run: /opt/sage/envs/sage/bin/sage code/claude_sib_curve3_lombardo.sage \
       > results/claude_sib_curve3_lombardo.log
"""

from sage.all import HyperellipticCurve, PolynomialRing, QQ

R = PolynomialRing(QQ, "x")
x = R.gen()

f = (
    3703062294195264 * x**6
    - 360079374491052216 * x**5
    + 8901721379573296848 * x**4
    - 5397945250386334945 * x**3
    - 86737535708373850908 * x**2
    + 36346694984390901540 * x
    + 43035470132681030400
)
h = x**2 + x

g = 4 * f + h**2
assert all(c.is_integer() for c in g.coefficients())

C = HyperellipticCurve(g)
J = C.jacobian()

alg_is_field = J.geometric_endomorphism_algebra_is_field(B=100)
ring_is_ZZ = J.geometric_endomorphism_ring_is_ZZ(B=100)

print("geometric_endomorphism_algebra_is_field(B=100) =", alg_is_field)
print("geometric_endomorphism_ring_is_ZZ(B=100) =", ring_is_ZZ)

if not alg_is_field:
    raise RuntimeError("FAIL: endomorphism algebra is NOT a field -> not geometrically simple")
if not ring_is_ZZ:
    raise RuntimeError("FAIL: endomorphism ring not certified ZZ")

print("LOMBARDO_CROSSCHECK_PASSED: End(Jac_Qbar) = ZZ, geometrically simple")
