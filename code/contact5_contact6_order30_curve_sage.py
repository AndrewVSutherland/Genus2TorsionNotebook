#!/usr/bin/env sage -python
"""Geometry probe for the simultaneous contact-5/contact-6 cover."""

from sage.all import QQ, PolynomialRing, Curve, ProjectiveSpace

R = PolynomialRing(QQ, ("u", "s"))
u, s = R.gens()

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

print("# contact5/contact6 order30 eliminated core")
print(f"degree_u={core.degree(u)} degree_s={core.degree(s)} total_degree={core.total_degree()} terms={len(core.monomials())}")
print(f"factorization={core.factor()}")

try:
    Caff = Curve(core)
    print(f"affine curve={Caff}")
    print(f"genus={Caff.genus()}")
except Exception as exc:  # noqa: BLE001
    print(f"genus_failed={type(exc).__name__}: {exc}")

try:
    Sring = PolynomialRing(QQ, "S")
    Sv = Sring.gen()
    pts = []
    for denu in range(1, 31):
        for numu in range(-30, 31):
            uu = QQ(numu) / QQ(denu)
            if uu in (0, 1, -1):
                continue
            pol_s = Sring(core.subs({u: uu, s: Sv}))
            if pol_s == 0:
                continue
            for rr, mult in pol_s.roots(QQ):
                pts.append((uu, rr))
    # Deduplicate while preserving order.
    seen = set()
    uniq = []
    for pt in pts:
        if pt not in seen:
            seen.add(pt)
            uniq.append(pt)
    print(f"rational fibers with height_u<=30: {len(uniq)}")
    for pt in uniq[:50]:
        print(f"point u={pt[0]} s={pt[1]}")
except Exception as exc:  # noqa: BLE001
    print(f"point_search_failed={type(exc).__name__}: {exc}")

try:
    P2 = ProjectiveSpace(QQ, 2, names=("U", "S", "Z"))
    U, S, Z = P2.coordinate_ring().gens()
    Rh = P2.coordinate_ring()
    hom = Rh(core.homogenize("Z"))
    Cproj = Curve(P2, hom)
    print(f"projective genus={Cproj.genus()}")
    sing = Cproj.singular_points()
    print(f"projective singular points count={len(sing)}")
    for P in sing[:30]:
        print(f"singular {P}")
except Exception as exc:  # noqa: BLE001
    print(f"projective_failed={type(exc).__name__}: {exc}")
