#!/usr/bin/env python3
"""Fast enumerator for the contact-7 first-halving surface.

This scans rational triples (s,u,z) and eps=+-1 using exact Fraction
arithmetic.  It writes candidate rows suitable for Magma verification:

    s u z eps a b r v w

The formulas match code/contact7_halving_surface_search.m.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd


def qstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


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


def coeffs_from_s_b(s: Fraction, b: Fraction, eps: int) -> tuple[Fraction, Fraction, Fraction, Fraction]:
    """Return c3,c2,c1,c0 for f(X+r)/X in the root chart."""
    S = Fraction(eps) * s
    den1 = S + 1
    if den1 == 0:
        raise ZeroDivisionError

    c3 = b * b - 5 * s * s - 2

    c2_num = (
        -2 * b * b * S**4
        - 4 * b * b * S**3
        + 4 * b * b * S
        + 2 * b * b
        + 2 * b * S**5
        + 4 * b * S**4
        + 6 * b * S**3
        + 8 * b * S**2
        + 10 * b * S
        + 5 * b
        + 10 * S**6
        + 20 * S**5
        + 18 * S**4
        + 16 * S**3
        + 11 * S**2
        + 6 * S
        + 3
    )
    c2 = c2_num / den1**2

    n1 = (
        -4 * b * b * S**8
        - 16 * b * b * S**7
        - 16 * b * b * S**6
        + 16 * b * b * S**5
        + 40 * b * b * S**4
        + 16 * b * b * S**3
        - 16 * b * b * S**2
        - 16 * b * b * S
        - 4 * b * b
        + 16 * b * S**9
        + 64 * b * S**8
        + 112 * b * S**7
        + 128 * b * S**6
        + 128 * b * S**5
        + 100 * b * S**4
        + 16 * b * S**3
        - 56 * b * S**2
        - 48 * b * S
        - 12 * b
        + 36 * S**10
        + 144 * S**9
        + 248 * S**8
        + 272 * S**7
        + 224 * S**6
        + 140 * S**5
        + 56 * S**4
        - 12 * S**3
        - 48 * S**2
        - 36 * S
        - 9
    )
    c1 = -n1 / (4 * den1**4)

    a0 = 3 * S**4 + 9 * S**3 + 11 * S**2 + 9 * S + 3
    c0 = S**7 * a0 / den1**3 + 2 * b * S**7
    return c3, c2, c1, c0


def b_from_s_z(s: Fraction, z: Fraction, eps: int) -> Fraction | None:
    S = Fraction(eps) * s
    if S == 0 or S == -1:
        return None
    den1 = S + 1
    a0 = 3 * S**4 + 9 * S**3 + 11 * S**2 + 9 * S + 3
    c00 = S**7 * a0 / den1**3
    lam = 2 * S**7
    if lam == 0:
        return None
    return (z * z - c00) / lam


def root_a_r(s: Fraction, b: Fraction, eps: int) -> tuple[Fraction, Fraction] | None:
    r = 1 - s * s
    if r == 0:
        return None
    a = (Fraction(eps) * s**7 - 1 + Fraction(7, 2) * r - b * r**3) / r**2
    return a, r


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--s-height", type=int, default=8)
    ap.add_argument("--u-height", type=int, default=8)
    ap.add_argument("--z-height", type=int, default=20)
    ap.add_argument("--out", default="")
    ap.add_argument("--progress", type=int, default=500000)
    args = ap.parse_args()

    s_vals = params(args.s_height)
    u_vals = params(args.u_height)
    z_vals = params(args.z_height)

    checked = 0
    hits = 0
    lines: list[str] = []

    for s in s_vals:
        if s == 0 or s * s == 1:
            continue
        for u in u_vals:
            if u == 0:
                continue
            for z in z_vals:
                for eps in (-1, 1):
                    checked += 1
                    if args.progress and checked % args.progress == 0:
                        print(f"progress checked={checked} hits={hits}", flush=True)
                    b = b_from_s_z(s, z, eps)
                    if b is None:
                        continue
                    c3, c2, c1, c0 = coeffs_from_s_b(s, b, eps)
                    if c0 != z * z:
                        continue
                    v = (u * u - c3) / 2
                    w = (v * v + 2 * z - c2) / (2 * u)
                    if c1 != w * w - 2 * v * z:
                        continue
                    ar = root_a_r(s, b, eps)
                    if ar is None:
                        continue
                    a, r = ar
                    hits += 1
                    lines.append(
                        " ".join(
                            qstr(q)
                            for q in (s, u, z, Fraction(eps), a, b, r, v, w)
                        )
                    )

    header = (
        f"# contact7 halving surface candidates\n"
        f"# s_height={args.s_height} u_height={args.u_height} z_height={args.z_height}\n"
        f"# checked={checked} hits={hits}\n"
        "# columns: s u z eps a b r v w\n"
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
