#!/usr/bin/env python3
"""
claude_ov_2216_hits2uv.py  -- lane OV/2216

Turn the C sweeper's "HIT <surface> <p> <q> <a> <b>" lines into exact
(u,v) parameters of the odd M_1(8,2,2) chart, re-verify EVERY square
condition in exact rational arithmetic (the C stage is only a modular
necessary condition), drop degenerate/singular parameters, and reduce
modulo the S_3 symmetry that permutes the triple {u, v, w=-(u+v+1)}.

Conditions re-verified here (see code/claude_ov_2216_symbolic.m):
  surface 0 (untwisted):  u(u-1)v(v-1),  u(u-1)(u+v+1)(u+v+2),  uv(u+v+1)
  surface 1 (twisted):   -u(u-1)v(2u+v+1), u(u-1)(u-v)(u+v+1),  uv(u+v+1)

Usage:  python3 claude_ov_2216_hits2uv.py <hitfile> [<hitfile> ...]
Prints one "u v" per surviving S_3 class, plus a summary.
"""
import sys
from fractions import Fraction
from math import isqrt, gcd


def is_sq(x: Fraction) -> bool:
    if x < 0:
        return False
    if x == 0:
        return False          # degenerate, treated separately
    n, d = x.numerator, x.denominator
    rn, rd = isqrt(n), isqrt(d)
    return rn * rn == n and rd * rd == d


def main():
    files = sys.argv[1:]
    if not files:
        print("usage: claude_ov_2216_hits2uv.py <hitfile> ...", file=sys.stderr)
        return 1

    seen = {}
    nlines = ndegen = nfail = 0
    for fn in files:
        for line in open(fn):
            if not line.startswith("HIT"):
                continue
            nlines += 1
            _, s, p, q, a, b = line.split()
            s, p, q, a, b = int(s), int(p), int(q), int(a), int(b)
            u = Fraction(p, q)
            A = u * (u - 1)
            k2 = Fraction(a * a, b * b)
            if s == 0:
                den = A - k2
                if den == 0:
                    ndegen += 1
                    continue
                v = A / den
            else:
                den = k2 + A
                if den == 0:
                    ndegen += 1
                    continue
                v = -A * (2 * u + 1) / den
            w = -(u + v + 1)
            # degeneracies: distinct sextic roots {1,u,v,w}, none equal
            vals = [Fraction(1), u, v, w]
            if len(set(vals)) != 4 or u == 0 or v == 0 or w == 0:
                ndegen += 1
                continue
            # exact re-verification
            c3 = u * v * (u + v + 1)
            if s == 0:
                c1 = u * (u - 1) * v * (v - 1)
                c2 = u * (u - 1) * (u + v + 1) * (u + v + 2)
            else:
                c1 = -u * (u - 1) * v * (2 * u + v + 1)
                c2 = u * (u - 1) * (u - v) * (u + v + 1)
            if not (is_sq(c1) and is_sq(c2) and is_sq(c3)):
                nfail += 1
                continue
            key = (s, tuple(sorted([u, v, w])))
            seen.setdefault(key, (u, v))

    print(f"# hit lines read      : {nlines}")
    print(f"# degenerate dropped  : {ndegen}")
    print(f"# failed exact recheck: {nfail}")
    print(f"# distinct S_3 classes: {len(seen)}")
    for (s, trip), (u, v) in sorted(seen.items(), key=lambda t: str(t[0])):
        print(f"{u} {v}    # surface {s}, triple {{{trip[0]}, {trip[1]}, {trip[2]}}}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
