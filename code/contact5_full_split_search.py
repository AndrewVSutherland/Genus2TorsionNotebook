#!/usr/bin/env python3
"""Search full rational 2-splitting in the contact-5/order-20 extra-2 locus.

On the double-linear locus of the scaled residual quartic, write p=r^2.  The
condition that two roots z,w are rational is

    Y^2 = (r+1)(r^2+2r+2)(r^3-r^2-4r+2).

The remaining quadratic factor splits exactly when

    W^2 = (r+2)(r^2+2r+2)(r^3-4r^2-2r+4).

For r=a/b in lowest terms both right hand sides have denominator b^6, so the
search is an integer square test on the two degree-6 numerators.
"""

from __future__ import annotations

import argparse
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path


def qstr(q: Fraction) -> str:
    return str(q.numerator) if q.denominator == 1 else f"{q.numerator}/{q.denominator}"


def is_square_int(n: int) -> int | None:
    if n < 0:
        return None
    r = isqrt(n)
    return r if r * r == n else None


def nums_for_ab(a: int, b: int) -> tuple[int, int]:
    q = a * a + 2 * a * b + 2 * b * b
    f1 = (a + b) * q * (a**3 - a * a * b - 4 * a * b * b + 2 * b**3)
    f2 = (a + 2 * b) * q * (a**3 - 4 * a * a * b - 2 * a * b * b + 4 * b**3)
    return f1, f2


def t_from_z(z: Fraction) -> Fraction | None:
    den = z**4 + 4 * z**3 + 8 * z * z + 8 * z + 4
    if den == 0:
        return None
    return -(z**4 + 4 * z + 4) / den


def roots_from_r_y(r: Fraction, y: Fraction) -> tuple[Fraction, Fraction] | None:
    if r == -2:
        return None
    s = -(r**3 + r * r + 2) / (r + 2)
    delta = y / (r + 2)
    return (s + delta) / 2, (s - delta) / 2


def residues_allowed(primes: list[int]) -> dict[int, set[int]]:
    out: dict[int, set[int]] = {}
    for p in primes:
        squares = {(i * i) % p for i in range(p)}
        allowed: set[int] = set()
        for r in range(p):
            q = (r * r + 2 * r + 2) % p
            f1 = ((r + 1) * q * (r**3 - r * r - 4 * r + 2)) % p
            f2 = ((r + 2) * q * (r**3 - 4 * r * r - 2 * r + 4)) % p
            if f1 in squares and f2 in squares:
                allowed.add(r)
        out[p] = allowed
    return out


def passes_residue_filter(a: int, b: int, allowed: dict[int, set[int]]) -> bool:
    for p, residues in allowed.items():
        bm = b % p
        if bm == 0:
            continue
        r = (a % p) * pow(bm, -1, p) % p
        if r not in residues:
            return False
    return True


def classify_hit(r: Fraction, y: Fraction) -> tuple[str, Fraction | None, Fraction | None, Fraction | None]:
    roots = roots_from_r_y(r, y)
    if roots is None:
        return "boundary", None, None, None
    z, w = roots
    if z == w:
        return "collision", z, w, None
    t1 = t_from_z(z)
    t2 = t_from_z(w)
    if t1 is None or t2 is None or t1 != t2:
        return "bad_t", z, w, None
    if t1 in {Fraction(-1), Fraction(-3)}:
        return "singular", z, w, t1
    return "candidate", z, w, t1


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=int, default=1000)
    ap.add_argument("--prime-bound", type=int, default=97)
    ap.add_argument("--out", default="")
    ap.add_argument("--progress", type=int, default=500000)
    args = ap.parse_args()

    small_primes = [
        p
        for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]
        if p <= args.prime_bound
    ]
    allowed = residues_allowed(small_primes)

    lines: list[str] = []
    checked = 0
    residue_survivors = 0
    f1_squares = 0
    both_squares = 0
    candidates = 0
    boundary = 0

    for b in range(1, args.height + 1):
        for a in range(-args.height, args.height + 1):
            if gcd(a, b) != 1:
                continue
            checked += 1
            if args.progress and checked % args.progress == 0:
                print(
                    "progress",
                    f"checked={checked}",
                    f"residue_survivors={residue_survivors}",
                    f"f1_squares={f1_squares}",
                    f"both_squares={both_squares}",
                    f"candidates={candidates}",
                    flush=True,
                )
            if not passes_residue_filter(a, b, allowed):
                continue
            residue_survivors += 1
            f1, f2 = nums_for_ab(a, b)
            y_num = is_square_int(f1)
            if y_num is None:
                continue
            f1_squares += 1
            w_num = is_square_int(f2)
            if w_num is None:
                continue
            both_squares += 1

            r = Fraction(a, b)
            y = Fraction(y_num, b**3)
            w = Fraction(w_num, b**3)
            status, z, zz, t = classify_hit(r, y)
            if status == "candidate":
                candidates += 1
            else:
                boundary += 1
            rows = [(y, w, status, z, zz, t)]
            if y_num != 0:
                status2, z2, zz2, t2 = classify_hit(r, -y)
                rows.append((-y, w, status2, z2, zz2, t2))
            for yy, ww, st, z0, z1, t0 in rows:
                lines.append(
                    " ".join(
                        [
                            qstr(r),
                            qstr(yy),
                            qstr(ww),
                            st,
                            qstr(z0) if z0 is not None else "NA",
                            qstr(z1) if z1 is not None else "NA",
                            qstr(t0) if t0 is not None else "NA",
                        ]
                    )
                )

    header = (
        "# contact5 full residual split search\n"
        f"# height={args.height} prime_bound={args.prime_bound} primes={len(small_primes)}\n"
        f"# checked={checked} residue_survivors={residue_survivors} "
        f"f1_squares={f1_squares} both_squares={both_squares} "
        f"candidates={candidates} boundary_or_degenerate={boundary}\n"
        "# columns: r Y W status z w t\n"
    )
    text = header + "\n".join(lines) + ("\n" if lines else "")
    if args.out:
        Path(args.out).write_text(text, encoding="ascii")
    else:
        print(text, end="")
    print(
        "DONE",
        f"checked={checked}",
        f"residue_survivors={residue_survivors}",
        f"f1_squares={f1_squares}",
        f"both_squares={both_squares}",
        f"candidates={candidates}",
        f"boundary_or_degenerate={boundary}",
    )


if __name__ == "__main__":
    main()
