#!/usr/bin/env python3
"""CRT reconstruction on the locally strongest full-quadratic-B slice c=2.

This imports the independent finite-field Hermite/M(12) matcher in
``m12_general5_fullquad_rootpair.py``.  For each requested prime it records
all smooth open matches with c=2, combines one local point at each prime, and
tries bounded rational reconstruction of every coordinate.  A reconstructed
tuple is accepted only after an exact Q-polynomial identity

    S_(s,t,a)(Z) = tau * F_(b,w)(d*Z-2)

with nonzero square tau.  Thus CRT coincidences cannot become false hits.

The default primes 7,11,13,17,19,23 give only 240 local combinations and use
negligible memory.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import math
import sys
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location(
    "fullquad", HERE / "m12_general5_fullquad_rootpair.py"
)
assert SPEC and SPEC.loader
FQ = importlib.util.module_from_spec(SPEC)
sys.modules["fullquad"] = FQ
SPEC.loader.exec_module(FQ)


def crt_pair(r: int, m: int, s: int, p: int) -> tuple[int, int]:
    z = ((s - r) * pow(m, -1, p)) % p
    return r + m * z, m * p


def crt(values: tuple[int, ...], primes: tuple[int, ...]) -> tuple[int, int]:
    r, m = 0, 1
    for s, p in zip(values, primes):
        r, m = crt_pair(r, m, s, p)
    return r, m


def best_reconstruction(r: int, modulus: int, bound: int) -> Fraction | None:
    """Lowest max(|numerator|,denominator) reduced lift within ``bound``."""
    best: tuple[int, int, int] | None = None
    for d in range(1, bound + 1):
        if math.gcd(d, modulus) != 1:
            continue
        n = (r * d) % modulus
        if n > modulus // 2:
            n -= modulus
        if abs(n) > bound or math.gcd(abs(n), d) != 1:
            continue
        row = (max(abs(n), d), abs(n) + d, d)
        if best is None or row < best:
            best = row
            ans = Fraction(n, d)
    return ans if best is not None else None


def trim(f):
    f = list(f)
    while len(f) > 1 and not f[-1]:
        f.pop()
    return f


def padd(f, g):
    n = max(len(f), len(g))
    return trim(
        [(f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)
         for i in range(n)]
    )


def pscale(f, a):
    return trim([a * x for x in f])


def pmul(f, g):
    h = [Fraction(0)] * (len(f) + len(g) - 1)
    for i, x in enumerate(f):
        for j, y in enumerate(g):
            h[i + j] += x * y
    return trim(h)


def ppow(f, n):
    ans, base = [Fraction(1)], list(f)
    while n:
        if n & 1:
            ans = pmul(ans, base)
        base = pmul(base, base)
        n >>= 1
    return ans


def pdiv_exact(f, g):
    r = trim(f)
    q = [Fraction(0)] * max(1, len(r) - len(g) + 1)
    while len(r) >= len(g) and r != [0]:
        z = r[-1] / g[-1]
        k = len(r) - len(g)
        q[k] = z
        for j, y in enumerate(g):
            r[j + k] -= z * y
        r = trim(r)
    if any(r):
        raise ArithmeticError(f"nonzero remainder {r}")
    return trim(q)


def hermite_exact(s, t, a):
    e = t * t - s * s - 1
    q = [s * s, e, Fraction(1)]
    A = [
        s**5,
        Fraction(5, 2) * s**3 * e,
        a + 2*s**5 - 5*s**3*t**2 + 5*s**3
        + Fraction(5, 2)*s**2*t**3 + Fraction(1, 2)*t**5
        - Fraction(5, 2)*t**3 + 2,
        -2*a - Fraction(1, 2)*s**5 + Fraction(5, 2)*s**3*t**2
        - Fraction(5, 2)*s**3 - Fraction(5, 2)*s**2*t**3
        + Fraction(1, 2)*t**5 + Fraction(5, 2)*t**3 - 3,
        a,
        Fraction(1),
    ]
    S = pdiv_exact(padd(pmul(A, A), pscale(ppow(q, 5), -1)), [0, 0, 1, -2, 1])
    return S, A, q


def compact_exact(b, w, d):
    # x=d*Z-2 (the fixed c=2 slice).
    x = [Fraction(-2), d]
    L = padd([b], pscale(x, 2*b - 1))
    H = padd(x, pscale(padd([1], pscale(x, b)), w))
    inner = padd(pmul(L, pmul(H, H)),
                 pscale(pmul(pmul(padd([1], x), padd([1], x)),
                             padd(pscale(L, w), pscale(pmul(x, x), -1))),
                        4*b))
    return pmul(L, inner)


def is_fraction_square(x: Fraction) -> bool:
    if x <= 0:
        return False
    return math.isqrt(x.numerator) ** 2 == x.numerator and \
        math.isqrt(x.denominator) ** 2 == x.denominator


def exact_match(vals: tuple[Fraction, ...]):
    b, w, d, s, t, a = vals
    if not b or b == 1 or 2*b == 1 or not w or not d or not s or not t:
        return None
    S, A, q = hermite_exact(s, t, a)
    G = compact_exact(b, w, d)
    S += [Fraction(0)] * (6 - len(S))
    G += [Fraction(0)] * (6 - len(G))
    if not G[5]:
        return None
    tau = S[5] / G[5]
    if not tau or any(S[i] != tau * G[i] for i in range(6)):
        return None
    if not is_fraction_square(tau):
        return None
    discq = q[1] ** 2 - 4*q[0]
    if not discq:
        return None
    return tau, A, q, S, G


def local_points(p: int):
    table, _ = FQ.build_S_hash(p)
    out = []
    c = 2 % p
    for b in range(p):
        for w in range(p):
            if not FQ.base_open(b, w, p):
                continue
            if not FQ.squarefree(FQ.compact_F_transformed(b, w, 0, 1, p), p):
                continue
            for d in range(1, p):
                G = FQ.compact_F_transformed(b, w, c, d, p)
                key, gpivot = FQ.projective_key(G, p)
                for rec in table.get(key, []):
                    ans = FQ.check_match(b, w, c, d, rec, G, gpivot, p)
                    if ans is not None:
                        # Keep b,w,d,s,t,a; c is the fixed constant 2.
                        v = ans[0]
                        out.append((v[0], v[1], v[3], v[4], v[5], v[6]))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", default="7,11,13,17,19,23")
    ap.add_argument("--bound", type=int, default=2000)
    ap.add_argument("--show", type=int, default=20,
                    help="number of lowest-height reconstructed rows to print")
    args = ap.parse_args()
    primes = tuple(int(x) for x in args.primes.split(",") if x)
    rows = []
    for p in primes:
        pts = local_points(p)
        print(f"LOCAL p={p} points={len(pts)} data={pts}")
        if not pts:
            print("C2_SLICE_LOCALLY_EMPTY")
            return
        rows.append(pts)
    modulus = math.prod(primes)
    ranked = []
    hits = []
    cache = {}
    for choice in itertools.product(*rows):
        vals = []
        for j in range(6):
            r, M = crt(tuple(x[j] for x in choice), primes)
            key = (r, M, args.bound)
            if key not in cache:
                cache[key] = best_reconstruction(r, M, args.bound)
            vals.append(cache[key])
        if any(x is None for x in vals):
            continue
        qvals = tuple(vals)
        score = max(max(abs(x.numerator), x.denominator) for x in qvals)
        ranked.append((score, qvals))
        ans = exact_match(qvals)
        if ans is not None:
            hits.append((qvals, ans[0]))
    ranked.sort(key=lambda z: (z[0], z[1]))
    print(f"CRT_SUMMARY primes={primes} modulus={modulus} "
          f"combinations={math.prod(len(x) for x in rows)} "
          f"fully_reconstructed={len(ranked)} exact_hits={len(hits)}")
    for score, vals in ranked[:args.show]:
        print(f"RECON score={score} b,w,d,s,t,a={vals}")
    for vals, tau in hits:
        print(f"EXACT_HIT b,w,d,s,t,a={vals} tau={tau}")


if __name__ == "__main__":
    main()
