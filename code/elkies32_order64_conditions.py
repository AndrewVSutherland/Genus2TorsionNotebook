#!/usr/bin/env sage -python
"""
Family-level halving conditions for the reconstructed Elkies N=32 family.

This derives conditions for D=(0,1)-infinity to be divisible by 2.  On the
reconstructed N=32 family this would give a point of order 64.
"""

from __future__ import annotations

from pathlib import Path
from sage.all import GF, PolynomialRing, QQ

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "data" / "elkies32_order64_conditions.txt"


def main():
    lines = []
    lines.append("# Elkies N=32 family: order-64 halving conditions")
    lines.append("")
    lines.append("For F(x)=sum F_i x^i and P=(0,1), D=P-infinity is divisible by 2 iff")
    lines.append("there are m,n and k=+/-1 such that, with")
    lines.append("  alpha = (F4 - m^2)/(2*F5),")
    lines.append("  beta  = (F3 - 2*m*n - F5*alpha^2)/(2*F5),")
    lines.append("the two equations")
    lines.append("  K2 = F2 - n^2 - 2*m*k - 2*F5*alpha*beta = 0,")
    lines.append("  K1 = F1 - 2*n*k - F5*beta^2 = 0")
    lines.append("hold.  Equivalently F(x)-(m*x^2+n*x+k)^2 = F5*x*(x^2+alpha*x+beta)^2.")
    lines.append("")

    R = PolynomialRing(QQ, ("z", "p", "m", "n"))
    z, p, m, n = R.gens()
    K = R.fraction_field()
    zk, pk, mk, nk = map(K, (z, p, m, n))

    a = zk * pk
    c = (-pk*zk**4 + 4*pk*zk**3 - 8*pk*zk**2 + 8*zk**2 + 8*pk*zk - 16)/(4*zk**2)
    b = pk*(zk - 1) + c - 1

    f1 = 2*c
    f2 = c**2 + 2*b
    f3 = 2*a + 2*b*c
    f4 = b**2 + 2*a*c
    f5 = a*(2*b - a)

    lines.append("After the first N=32 rank condition and z=a/p:")
    lines.append("  a = z*p")
    lines.append(f"  c = {c}")
    lines.append(f"  b = {b}")
    lines.append("  F1=2*c, F2=c^2+2*b, F3=2*a+2*b*c, F4=b^2+2*a*c, F5=a*(2*b-a)")
    lines.append("")

    for sign in [-1, 1]:
        alpha = (f4 - mk**2)/(2*f5)
        beta = (f3 - 2*mk*nk - f5*alpha**2)/(2*f5)
        K2 = f2 - nk**2 - 2*mk*sign - 2*f5*alpha*beta
        K1 = f1 - 2*nk*sign - f5*beta**2
        K2n = R(K2.numerator())
        K1n = R(K1.numerator())
        lines.append(f"Halving sign k={sign} after clearing denominators:")
        lines.append(f"  K2 numerator degrees in (z,p,m,n): {[K2n.degree(v) for v in (z,p,m,n)]}; terms={len(K2n.dict())}")
        lines.append(f"  K1 numerator degrees in (z,p,m,n): {[K1n.degree(v) for v in (z,p,m,n)]}; terms={len(K1n.dict())}")
    lines.append("")

    lines.append("On the reconstructed N=32 component, (z,r) additionally satisfies")
    lines.append("  C32(z,r)=0, the bidegree (4,4) genus-0 equation in data/elkies32_reconstruct_conditions.txt,")
    lines.append("and p is the rational common-root formula recorded there.  Thus the order-64 cover is")
    lines.append("  C32(z,r)=0, p=4*Pnum(z,r)/Pden(z,r), K1(z,p,m,n)=K2(z,p,m,n)=0")
    lines.append("for one of the two signs k=+/-1, away from the boundary denominators z*p*F5*Pden=0.")
    lines.append("")

    def Cval(zv, rv):
        return (3*zv**4*rv**4 + 9*zv**4*rv**3 - 16*zv**3*rv**4 + 10*zv**4*rv**2
                - 56*zv**3*rv**3 + 32*zv**2*rv**4 + 5*zv**4*rv - 72*zv**3*rv**2
                + 144*zv**2*rv**3 - 64*zv*rv**4 + zv**4 - 40*zv**3*rv
                + 208*zv**2*rv**2 - 224*zv*rv**3 + 80*rv**4 - 8*zv**3
                + 120*zv**2*rv - 288*zv*rv**2 + 160*rv**3 + 24*zv**2
                - 160*zv*rv + 160*rv**2 - 32*zv + 80*rv + 16)

    def Pnum(zv, rv):
        return ((zv - 2)*(rv + 1)**2*(zv**2*rv + zv**2 - 4*rv)
                *(zv**4*rv**2 + 2*zv**4*rv - 10*zv**3*rv**2 + zv**4
                  - 16*zv**3*rv + 28*zv**2*rv**2 - 6*zv**3 + 44*zv**2*rv
                  - 16*zv*rv**2 + 16*zv**2 - 56*zv*rv + 16*rv**2
                  - 24*zv + 32*rv + 16))

    def Pden(zv, rv):
        return (rv*zv*(zv**8*rv**4 + 4*zv**8*rv**3 - 24*zv**7*rv**4
                + 48*zv**6*rv**5 + 6*zv**8*rv**2 - 80*zv**7*rv**3
                + 284*zv**6*rv**4 - 256*zv**5*rv**5 + 4*zv**8*rv
                - 100*zv**7*rv**2 + 632*zv**6*rv**3 - 1200*zv**5*rv**4
                + 320*zv**4*rv**5 + zv**8 - 56*zv**7*rv + 672*zv**6*rv**2
                - 2432*zv**5*rv**3 + 2224*zv**4*rv**4 - 12*zv**7
                + 344*zv**6*rv - 2496*zv**5*rv**2 + 5312*zv**4*rv**3
                - 2304*zv**3*rv**4 + 68*zv**6 - 1248*zv**5*rv
                + 5808*zv**4*rv**2 - 7232*zv**3*rv**3 + 1664*zv**2*rv**4
                - 240*zv**5 + 2976*zv**4*rv - 8832*zv**3*rv**2
                + 6272*zv**2*rv**3 - 768*zv*rv**4 + 576*zv**4
                - 4800*zv**3*rv + 8640*zv**2*rv**2 - 3328*zv*rv**3
                + 256*rv**4 - 960*zv**3 + 5120*zv**2*rv - 5120*zv*rv**2
                + 1024*rv**3 + 1088*zv**2 - 3328*zv*rv + 1536*rv**2
                - 768*zv + 1024*rv + 256))

    def has_half_mod(F, zv, rv, pv):
        if zv == 0 or pv == 0:
            return False
        av = zv*pv
        cv = (-pv*zv**4 + 4*pv*zv**3 - 8*pv*zv**2 + 8*zv**2 + 8*pv*zv - 16)/(4*zv**2)
        bv = pv*(zv - 1) + cv - 1
        F1 = 2*cv
        F2 = cv**2 + 2*bv
        F3 = 2*av + 2*bv*cv
        F4 = bv**2 + 2*av*cv
        F5 = av*(2*bv - av)
        if F5 == 0:
            return False
        for sgn in (F(-1), F(1)):
            for mv in F:
                for nv in F:
                    al = (F4 - mv**2)/(2*F5)
                    be = (F3 - 2*mv*nv - F5*al**2)/(2*F5)
                    if F2 - nv**2 - 2*mv*sgn - 2*F5*al*be == 0 and F1 - 2*nv*sgn - F5*be**2 == 0:
                        return True
        return False

    lines.append("Affine finite-field sanity check on the derived cover:")
    lines.append("  q: #C32_affine  #usable_p_formula  #usable_with_halving")
    for q in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73]:
        Fq = GF(q)
        base = usable = half = 0
        for zv in Fq:
            for rv in Fq:
                if Cval(zv, rv) != 0:
                    continue
                base += 1
                den = Pden(zv, rv)
                if den == 0 or zv == 0:
                    continue
                pv = 4*Pnum(zv, rv)/den
                if pv == 0:
                    continue
                usable += 1
                if has_half_mod(Fq, zv, rv, pv):
                    half += 1
        lines.append(f"  {q}: {base:3d} {usable:3d} {half:3d}")

    OUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
