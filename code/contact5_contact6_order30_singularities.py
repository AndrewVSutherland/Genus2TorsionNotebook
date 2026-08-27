#!/usr/bin/env sage -python
"""Projective singularity probe for the order-30 eliminated core."""

from sage.all import QQ, PolynomialRing, ProjectiveSpace, Curve

R2 = PolynomialRing(QQ, ("u", "s"))
u, s = R2.gens()
core = (
    u**10 + 20*u**9 + 6*u**8*s + 223*u**8 + 30*u**7*s
    - 21*u**6*s**2 - 2*u**5*s**3 + 1380*u**7
    - 438*u**6*s - 60*u**5*s**2 + 34*u**4*s**3
    - 6*u**3*s**4 + 4005*u**6 - 3525*u**5*s
    + 1557*u**4*s**2 - 488*u**3*s**3 + 105*u**2*s**4
    - 15*u*s**5 + s**6 + 2796*u**5 - 2256*u**4*s
    + 780*u**3*s**2 - 150*u**2*s**3 + 12*u*s**4
    + 767*u**4 - 420*u**3*s + 75*u**2*s**2 - 4*u*s**3
    + 70*u**3 - 12*u**2*s - u**2
)

R3 = PolynomialRing(QQ, ("U", "S", "Z"))
U, S, Z = R3.gens()
hom = R3(0)
for (i, j), coeff in core.dict().items():
    hom += QQ(coeff) * U**i * S**j * Z**(10 - i - j)

P2 = ProjectiveSpace(QQ, 2, names=("U", "S", "Z"))
U2, S2, Z2 = P2.coordinate_ring().gens()
hom2 = P2.coordinate_ring()(hom(U2, S2, Z2))
C = Curve(hom2)
print(f"degree={C.degree()} genus={C.genus()}")
try:
    sings = C.singular_points()
    print(f"rational singular points={len(sings)}")
    for P in sings:
        print(P)
except Exception as exc:  # noqa: BLE001
    print(f"singular_points_failed={type(exc).__name__}: {exc}")

Fx = hom.derivative(U)
Fy = hom.derivative(S)
Fz = hom.derivative(Z)
for P in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, -1, 0), (1, 1, 0)]:
    vals = [g(*P) for g in (hom, Fx, Fy, Fz)]
    print(f"P={P} vals={vals}")

def multiplicity_at_affine(poly, x0, y0):
    T = PolynomialRing(QQ, ("X", "Y"))
    X, Y = T.gens()
    g = T(poly(X + x0, Y + y0))
    if g == 0:
        return -1
    return min(sum(exp) for exp, coeff in g.dict().items() if coeff)


for Paff in [(-1, -3), (0, 0), (QQ(1)/8, QQ(21)/32)]:
    print(f"affine point {Paff} multiplicity", multiplicity_at_affine(core, Paff[0], Paff[1]))

# Infinity charts.
RSZ = PolynomialRing(QQ, ("v", "z"))
v, z = RSZ.gens()
chart_U = RSZ(hom(1, v, z))
print("chart U=1 at (S,Z)=(0,0) multiplicity", multiplicity_at_affine(chart_U, 0, 0))
print("chart U=1 polynomial low terms", chart_U)

RUZ = PolynomialRing(QQ, ("w", "z"))
w, z2 = RUZ.gens()
chart_S = RUZ(hom(w, 1, z2))
print("chart S=1 at (U,Z)=(0,0) multiplicity", multiplicity_at_affine(chart_S, 0, 0))
