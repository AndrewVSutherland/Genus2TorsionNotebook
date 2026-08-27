#!/usr/bin/env python3
"""Enumerate simultaneous contact-7 and point-contact-5 candidates.

The contact-7 family is

    f = (h7^2 + (x-1)^7)/x^2,
    h7 = 1 - 7*x/2 + a*x^2 + b*x^3.

We impose a point-contact 5 condition

    q^2 - f = const*(x-r)^5,   q = c0 + c1*x + c2*x^2.

The top derivative equation is

    c2^2 - b^2 = 5*r - 7.

Write c2-b=d and c2+b=e, so r=(d*e+7)/5,
b=(e-d)/2, c2=(e+d)/2.  Then the next two derivative
equations solve for a and c0, leaving the first two equations as checks.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd


def params(height: int) -> list[Fraction]:
    out: list[Fraction] = []
    seen: set[Fraction] = set()
    for den in range(1, height + 1):
        for num in range(-height, height + 1):
            if gcd(num, den) != 1:
                continue
            q = Fraction(num, den)
            if q not in seen:
                seen.add(q)
                out.append(q)
    return out


def qstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def equations_ok(a: Fraction, b: Fraction, c0: Fraction, c1: Fraction, c2: Fraction, r: Fraction) -> bool:
    e0 = -(
        4*a*a*r*r + 8*a*b*r**3 - 28*a*r + 8*a
        + 4*b*b*r**4 - 28*b*r*r + 8*b*r
        - 4*c0*c0 - 8*c0*c1*r - 8*c0*c2*r*r
        - 4*c1*c1*r*r - 8*c1*c2*r**3 - 4*c2*c2*r**4
        + 4*r**5 - 28*r**4 + 84*r**3 - 140*r*r + 140*r - 35
    ) / 4
    if e0 != 0:
        return False

    e1 = (
        -2*a*a*r - 6*a*b*r*r + 7*a - 4*b*b*r**3 + 14*b*r - 2*b
        + 2*c0*c1 + 4*c0*c2*r + 2*c1*c1*r + 6*c1*c2*r*r
        + 4*c2*c2*r**3 - 5*r**4 + 28*r**3 - 63*r*r + 70*r - 35
    )
    return e1 == 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=int, default=8)
    ap.add_argument("--out", default="")
    ap.add_argument("--progress", type=int, default=250000)
    args = ap.parse_args()

    vals = params(args.height)
    checked = 0
    hits = 0
    lines: list[str] = []

    for d in vals:
        if d == 0:
            continue
        for e in vals:
            if e == 0:
                continue
            r = (d * e + 7) / 5
            b = (e - d) / 2
            c2 = (e + d) / 2
            if b == 0 or c2 == 0:
                continue
            for c1 in vals:
                checked += 1
                if args.progress and checked % args.progress == 0:
                    print(f"progress checked={checked} hits={hits}", flush=True)

                a_num = (
                    2*c1*c2 + 4*c2*c2*r - 4*b*b*r
                    - 10*r*r + 28*r - 21
                )
                a = a_num / (2*b)

                c0_num = (
                    a*a + 6*a*b*r + 6*b*b*r*r - 7*b - c1*c1
                    - 6*c1*c2*r - 6*c2*c2*r*r
                    + 10*r**3 - 42*r*r + 63*r - 35
                )
                c0 = c0_num / (2*c2)

                if not equations_ok(a, b, c0, c1, c2, r):
                    continue
                if r == 1:
                    continue
                hits += 1
                lines.append(" ".join(qstr(q) for q in (a, b, r, c0, c1, c2, d, e)))

    header = (
        "# simultaneous contact7/contact5 point candidates\n"
        f"# height={args.height} checked={checked} hits={hits}\n"
        "# columns: a b r c0 c1 c2 d e\n"
    )
    text = header + "\n".join(lines) + ("\n" if lines else "")
    if args.out:
        with open(args.out, "w", encoding="ascii") as fh:
            fh.write(text)
    else:
        print(text, end="")
    print(f"DONE checked={checked} hits={hits}")


if __name__ == "__main__":
    main()
