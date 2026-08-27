#!/usr/bin/env python3
"""Exact finite-support masks for the transverse HLP deformation slice.

The slice is

    F_t = F_HLP + t*(1+x).

For each finite t over F_p this script tests:

1. at least two distinct projective quadratic cubic-contact supports
       H^2-F_t = k*Q^3;
2. a projective quadratic factor q0 and a half support u satisfying
       (q0*L)^2-F_t = k*u^2*q0.

Degree-one/degree-zero affine representatives are included for projective
supports meeting infinity.  Singular curve residues are retained as boundary;
the projective parameter t=infinity is to be retained by the height sieve.
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path


def trim(a: list[int], p: int) -> list[int]:
    a = [x % p for x in a]
    while a and a[-1] == 0:
        a.pop()
    return a


def add(a: list[int], b: list[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) +
                 (b[i] if i < len(b) else 0) for i in range(n)], p)


def sub(a: list[int], b: list[int], p: int) -> list[int]:
    n = max(len(a), len(b))
    return trim([(a[i] if i < len(a) else 0) -
                 (b[i] if i < len(b) else 0) for i in range(n)], p)


def scale(a: list[int], c: int, p: int) -> list[int]:
    return trim([c * x for x in a], p)


def mul(a: list[int], b: list[int], p: int) -> list[int]:
    if not a or not b:
        return []
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % p
    return trim(out, p)


def power(a: list[int], n: int, p: int) -> list[int]:
    out = [1]
    for _ in range(n):
        out = mul(out, a, p)
    return out


def divrem(a: list[int], b: list[int], p: int) -> tuple[list[int], list[int]]:
    a, b = trim(a, p), trim(b, p)
    quo = [0] * max(1, len(a) - len(b) + 1)
    ib = pow(b[-1], -1, p)
    while a and len(a) >= len(b):
        j = len(a) - len(b)
        c = a[-1] * ib % p
        quo[j] = c
        for i, bi in enumerate(b):
            a[i + j] = (a[i + j] - c * bi) % p
        a = trim(a, p)
    return trim(quo, p), a


def gcd_poly(a: list[int], b: list[int], p: int) -> list[int]:
    while trim(b, p):
        a, b = b, divrem(a, b, p)[1]
    a = trim(a, p)
    return scale(a, pow(a[-1], -1, p), p) if a else []


def derivative(a: list[int], p: int) -> list[int]:
    return trim([i * a[i] for i in range(1, len(a))], p)


def pad7(a: list[int], p: int) -> list[int]:
    a = trim(a, p)
    return a + [0] * (7 - len(a))


def slice_polynomial(t: int, p: int) -> list[int]:
    return [
        (187392 + t) % p,
        t % p,
        -118767 % p,
        0,
        -118767 % p,
        0,
        187392 % p,
    ]


def good_sextic(f: list[int], p: int) -> bool:
    return len(trim(f, p)) - 1 == 6 and \
        len(gcd_poly(f, derivative(f, p), p)) == 1


def monic_contact_supports(f: list[int], p: int) -> set[tuple[int, int, int]]:
    """Supports (1,U,V) for Q=x^2+U*x+V."""
    f0, f1, f2, f3, f4, f5, f6 = pad7(f, p)
    supports: set[tuple[int, int, int]] = set()
    for U in range(p):
        for V in range(p):
            q3 = [
                V**3,
                3 * U * V**2,
                3 * (U**2 * V + V**2),
                U**3 + 6 * U * V,
                3 * (U**2 + V),
                3 * U,
                1,
            ]
            q3 = [x % p for x in q3]
            found = False
            for h3 in range(p):
                k = (h3 * h3 - f6) % p
                if h3:
                    inv2h3 = pow(2 * h3, -1, p)
                    h2 = (f5 + k * q3[5]) * inv2h3 % p
                    h1 = (f4 + k * q3[4] - h2 * h2) * inv2h3 % p
                    h0 = (f3 + k * q3[3] - 2 * h2 * h1) * inv2h3 % p
                    if (h1*h1 + 2*h2*h0 - f2 - k*q3[2]) % p:
                        continue
                    if (2*h1*h0 - f1 - k*q3[1]) % p:
                        continue
                    if (h0*h0 - f0 - k*q3[0]) % p:
                        continue
                    found = True
                    break

                # Exceptional h3=0 branch.
                if (-f5 - k*q3[5]) % p:
                    continue
                for h2 in range(p):
                    if (h2*h2 - f4 - k*q3[4]) % p:
                        continue
                    h1_values = range(p) if h2 == 0 else [
                        (f3 + k*q3[3]) * pow(2*h2, -1, p) % p
                    ]
                    for h1 in h1_values:
                        if h2:
                            h0_values = [
                                (f2 + k*q3[2] - h1*h1)
                                * pow(2*h2, -1, p) % p
                            ]
                        elif (h1*h1 - f2 - k*q3[2]) % p == 0:
                            h0_values = range(p)
                        else:
                            continue
                        for h0 in h0_values:
                            if (2*h1*h0 - f1 - k*q3[1]) % p:
                                continue
                            if (h0*h0 - f0 - k*q3[0]) % p:
                                continue
                            found = True
                            break
                        if found:
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                supports.add((1, U, V))
    return supports


def lower_degree_contact_supports(f: list[int], p: int) -> set[tuple[int, int, int]]:
    """Projective supports (0,1,V) and (0,0,1), i.e. supports at infinity."""
    f0, f1, f2, f3, f4, f5, f6 = pad7(f, p)
    supports: set[tuple[int, int, int]] = set()
    square_roots = [a for a in range(p) if a*a % p == f6]

    # Q=x+V, so Q^3 has coefficients V^3,3V^2,3V,1.
    for V in range(p):
        found = False
        for h3 in square_roots:
            h2 = f5 * pow(2*h3, -1, p) % p
            h1 = (f4 - h2*h2) * pow(2*h3, -1, p) % p
            for h0 in range(p):
                k = (2*h3*h0 + 2*h2*h1 - f3) % p
                if (h1*h1 + 2*h2*h0 - f2 - 3*k*V) % p:
                    continue
                if (2*h1*h0 - f1 - 3*k*V*V) % p:
                    continue
                if (h0*h0 - f0 - k*V**3) % p:
                    continue
                found = True
                break
            if found:
                break
        if found:
            supports.add((0, 1, V))

    # Q=1.  H^2-F must be constant.
    found = False
    for h3 in square_roots:
        h2 = f5 * pow(2*h3, -1, p) % p
        h1 = (f4 - h2*h2) * pow(2*h3, -1, p) % p
        h0 = (f3 - 2*h2*h1) * pow(2*h3, -1, p) % p
        if (h1*h1 + 2*h2*h0 - f2) % p:
            continue
        if (2*h1*h0 - f1) % p:
            continue
        found = True
        break
    if found:
        supports.add((0, 0, 1))
    return supports


def contact_supports(f: list[int], p: int) -> set[tuple[int, int, int]]:
    return monic_contact_supports(f, p) | lower_degree_contact_supports(f, p)


def scalar_quadratic_square(s: list[int], p: int):
    """Return (projective quadratic u,k) if s=k*u^2 with k nonzero."""
    s = pad7(s, p)[:5]
    s0, s1, s2, s3, s4 = s
    if s4:
        k = s4
        u1 = s3 * pow(2*k, -1, p) % p
        u0 = (s2 * pow(k, -1, p) - u1*u1) * pow(2, -1, p) % p
        if s1 == 2*k*u0*u1 % p and s0 == k*u0*u0 % p:
            return (1, u1, u0), k
        return None
    if s3:
        return None
    if s2:
        k = s2
        u0 = s1 * pow(2*k, -1, p) % p
        if s0 == k*u0*u0 % p:
            return (0, 1, u0), k
        return None
    if s1:
        return None
    if s0:
        return (0, 0, 1), s0
    return None


def halving_supports(f: list[int], p: int) -> set[tuple]:
    """Enumerate exact projective q0/u/L halving supports."""
    supports: set[tuple] = set()
    for c1 in range(p):
        for c0 in range(p):
            q0 = [c0, c1, 1]
            quotient, remainder = divrem(f, q0, p)
            if remainder:
                continue
            for l1 in range(p):
                for l0 in range(p):
                    line = [l0, l1]
                    s = sub(mul(q0, power(line, 2, p), p), quotient, p)
                    ans = scalar_quadratic_square(s, p)
                    if ans is None:
                        continue
                    u, k = ans
                    supports.add(((1, c1, c0), u, (l1, l0), k))
    return supports


def prime_mask(p: int):
    allowed = []
    bad = []
    contact_ok = []
    half_ok = []
    counts = {}
    for t in range(p):
        f = slice_polynomial(t, p)
        if not good_sextic(f, p):
            bad.append(t)
            allowed.append(t)
            continue
        contacts = contact_supports(f, p)
        halves = halving_supports(f, p)
        counts[t] = (len(contacts), len(halves))
        if len(contacts) >= 2:
            contact_ok.append(t)
        if halves:
            half_ok.append(t)
        if len(contacts) >= 2 and halves:
            allowed.append(t)
    return allowed, bad, contact_ok, half_ok, counts


def primes_up_to(n: int) -> list[int]:
    return [p for p in range(2, n + 1)
            if all(p % d for d in range(2, math.isqrt(p) + 1))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes")
    ap.add_argument("--prime-bound", type=int, default=43)
    ap.add_argument("--write-masks", type=Path)
    args = ap.parse_args()
    primes = ([int(x) for x in args.primes.split(",") if x]
              if args.primes else
              [p for p in primes_up_to(args.prime_bound) if p not in (2, 3)])
    lines = [
        "# p : simultaneous_allowed ; bad ; contact_ok ; halving_ok ; counts t=c,h",
        "# Projective t=infinity is always conservatively allowed.",
    ]
    for p in primes:
        allowed, bad, contacts, halves, counts = prime_mask(p)
        count_text = ",".join(f"{t}={c}/{h}" for t, (c, h) in sorted(counts.items()))
        lines.append(
            f"{p} : {','.join(map(str,allowed))} ; {','.join(map(str,bad))} ; "
            f"{','.join(map(str,contacts))} ; {','.join(map(str,halves))} ; {count_text}"
        )
        print(f"p={p} allowed={allowed} bad={bad} contact_ok={contacts} "
              f"halving_ok={halves} density={(len(allowed)+1)/(p+1):.6f}")
    if args.write_masks:
        args.write_masks.parent.mkdir(parents=True, exist_ok=True)
        args.write_masks.write_text("\n".join(lines) + "\n")
        print("wrote", args.write_masks)


if __name__ == "__main__":
    main()
