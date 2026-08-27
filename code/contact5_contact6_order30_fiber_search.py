#!/usr/bin/env sage -python
"""Search rational fibers on the genus-zero contact5/contact6 core."""

from __future__ import annotations

import argparse
from math import gcd

from sage.all import QQ, PolynomialRing


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=int, default=300)
    ap.add_argument("--out", default="data/contact5_contact6_order30_fibers_h300.txt")
    args = ap.parse_args()

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
    Sring = PolynomialRing(QQ, "S")
    Sv = Sring.gen()

    hits = []
    seen = set()
    checked = 0
    for den in range(1, args.height + 1):
        for num in range(-args.height, args.height + 1):
            if gcd(num, den) != 1:
                continue
            uu = QQ(num) / QQ(den)
            if uu in (0, 1, -1):
                continue
            checked += 1
            pol = Sring(core.subs({u: uu, s: Sv}))
            if pol == 0:
                continue
            for ss, _ in pol.roots(QQ):
                pt = (uu, ss)
                if pt in seen:
                    continue
                seen.add(pt)
                hits.append(pt)

    lines = [
        "# rational fibers on contact5/contact6 order30 core",
        f"# height={args.height} checked_u={checked} hits={len(hits)}",
        "# columns: u s",
    ]
    for uu, ss in hits:
        lines.append(f"{uu} {ss}")
    text = "\n".join(lines) + "\n"
    with open(args.out, "w", encoding="ascii") as fh:
        fh.write(text)
    print(text, end="")


if __name__ == "__main__":
    main()
