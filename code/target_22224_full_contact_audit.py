#!/usr/bin/env python3
"""Audit finite-field cubic-contact masks for the full A(2,2,2,8) cover.

For

    f = x(x+a^2)(x+b^2)(x+c^2)(x+d^2)

we enumerate identities

    h^2 - f = m^2 q^3,
    q = x^2 + U*x + v^2,  m = 1/L.

The formulas here are derived directly by coefficient comparison and every
generated coefficient key is checked by reconstructing the polynomial
identity.  This script exists separately from the earlier finite-mask code
because that code had a factor-of-two error in the x^4 comparison.
"""

from __future__ import annotations

import argparse
from collections import Counter


def add_poly(a, b, p):
    out = [0] * max(len(a), len(b))
    for i, z in enumerate(a):
        out[i] = (out[i] + z) % p
    for i, z in enumerate(b):
        out[i] = (out[i] + z) % p
    return out


def sub_poly(a, b, p):
    out = [0] * max(len(a), len(b))
    for i, z in enumerate(a):
        out[i] = (out[i] + z) % p
    for i, z in enumerate(b):
        out[i] = (out[i] - z) % p
    return out


def mul_poly(a, b, p):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] = (out[i + j] + x * y) % p
    return out


def square_set(p):
    return {x * x % p for x in range(1, p)}


def cover_radicands(t, p):
    a, b, c, d = t
    return (
        a * b * c * d % p,
        a * (a + b) * (a + c) * (a + d) % p,
        b * (b + a) * (b + c) * (b + d) % p,
        c * (c + a) * (c + b) * (c + d) % p,
    )


def projective_tuples(p):
    # Unique representative: the first nonzero coordinate is 1.
    for first in range(4):
        prefix = [0] * first + [1]
        tail = 3 - first
        if tail == 0:
            yield tuple(prefix)
        elif tail == 1:
            for x in range(p):
                yield tuple(prefix + [x])
        elif tail == 2:
            for x in range(p):
                for y in range(p):
                    yield tuple(prefix + [x, y])
        else:
            for x in range(p):
                for y in range(p):
                    for z in range(p):
                        yield tuple(prefix + [x, y, z])


def elementary_from_signed(t, p):
    A, B, C, D = (z * z % p for z in t)
    e1 = (A + B + C + D) % p
    e2 = (A * B + A * C + A * D + B * C + B * D + C * D) % p
    e3 = (A * B * C + A * B * D + A * C * D + B * C * D) % p
    e4 = A * B * C * D % p
    return e1, e2, e3, e4


def verify_identity(p, key, witness):
    e1, e2, e3, e4 = key
    L, U, v = witness
    inv = lambda z: pow(z % p, -1, p)
    M = L * L % p
    P = (4 * M * e1 + 12 * (U * U + v * v) - (M + 3 * U) ** 2) % p
    m = inv(L)
    A = (M + 3 * U) * inv(2 * L) % p
    B = P * inv(8 * L) % p
    C = v**3 * m % p
    h = [C, B, A, m]
    q = [v * v % p, U, 1]
    f = [0, e4, e3, e2, e1, 1]
    lhs = sub_poly(mul_poly(h, h, p), f, p)
    q3 = mul_poly(mul_poly(q, q, p), q, p)
    rhs = [(m * m * z) % p for z in q3]
    lhs += [0] * (7 - len(lhs))
    rhs += [0] * (7 - len(rhs))
    return all((x - y) % p == 0 for x, y in zip(lhs, rhs))


def contact_keys(p):
    """Return all coefficient keys and one exactly audited witness per key."""
    inv = lambda z: pow(z % p, -1, p)
    witnesses = {}
    discr_hist = Counter()
    for L in range(1, p):
        M = L * L % p
        i8M = inv(8 * M)
        i64M = inv(64 * M)
        i4M = inv(4 * M)
        for U in range(p):
            for v in range(p):
                for e1 in range(p):
                    P = (
                        4 * M * e1
                        + 12 * (U * U + v * v)
                        - (M + 3 * U) ** 2
                    ) % p
                    e2 = (
                        (M + 3 * U) * P
                        + 16 * v**3
                        - 8 * U**3
                        - 48 * U * v * v
                    ) * i8M % p
                    e3 = (
                        P * P
                        + 64 * (M + 3 * U) * v**3
                        - 192 * (U * U * v * v + v**4)
                    ) * i64M % p
                    e4 = (P * v**3 - 12 * U * v**4) * i4M % p
                    if e4 == 0:
                        continue
                    key = (e1, e2, e3, e4)
                    if key not in witnesses:
                        witness = (L, U, v)
                        if not verify_identity(p, key, witness):
                            raise AssertionError((p, key, witness))
                        witnesses[key] = witness
                        discr_hist[(U * U - 4 * v * v) % p == 0] += 1
    return witnesses, discr_hist


def audit_prime(p):
    squares = square_set(p)
    witnesses, discr_hist = contact_keys(p)
    counts = Counter()
    target = []
    target_keys = set()
    for t in projective_tuples(p):
        counts["projective"] += 1
        rad = cover_radicands(t, p)
        if not all(z in squares for z in rad):
            continue
        counts["full_cover"] += 1
        sq = [z * z % p for z in t]
        if 0 in sq or len(set(sq)) != 4:
            continue
        counts["smooth_full_cover"] += 1
        key = elementary_from_signed(t, p)
        if key not in witnesses:
            continue
        counts["target_presentations"] += 1
        target_keys.add(key)
        if len(target) < 12:
            target.append((t, key, witnesses[key]))
    counts["contact_keys"] = len(witnesses)
    counts["target_curvekeys"] = len(target_keys)
    return counts, discr_hist, target


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="5,7,11,13,17,19,23,29,31")
    args = parser.parse_args()
    for p in map(int, args.primes.split(",")):
        counts, discr_hist, samples = audit_prime(p)
        print("AUDIT", p, dict(counts), "contact_qsquare_hist", dict(discr_hist))
        for row in samples:
            print("SAMPLE", p, row)


if __name__ == "__main__":
    main()
