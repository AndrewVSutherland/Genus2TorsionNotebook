#!/usr/bin/env python3
"""Independent finite-field check of the reduced b2=0 sign quotient.

This uses the three generic equations K2,K3,K4 and the small exceptional
G(0)=0 system from m12_general5_b2zero_rootquotient.m.  It counts quotient
points which lift to rational b1 (equivalently tau is a nonzero square).
"""

from argparse import ArgumentParser
from math import comb


def trim(a):
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def add(a, b, p):
    r = [0] * max(len(a), len(b))
    for i in range(len(r)):
        r[i] = ((a[i] if i < len(a) else 0)
                + (b[i] if i < len(b) else 0)) % p
    return trim(r)


def scale(a, s, p):
    return trim([(s * x) % p for x in a])


def mul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            r[i + j] = (r[i + j] + x * y) % p
    return trim(r)


def mod_poly(a, b, p):
    a = trim(a[:])
    b = trim(b[:])
    ib = pow(b[-1], -1, p)
    while len(a) >= len(b) and a != [0]:
        s = a[-1] * ib % p
        k = len(a) - len(b)
        for i, x in enumerate(b):
            a[i + k] = (a[i + k] - s * x) % p
        trim(a)
    return a


def coprime(a, b, p):
    a, b = trim(a[:]), trim(b[:])
    while b != [0]:
        a, b = b, mod_poly(a, b, p)
    return len(a) == 1 and a[0] != 0


def compact_f(b, w, p):
    L = [b, (2 * b - 1) % p]
    H = [w, (1 + w * b) % p]
    one_x_sq = [1, 2, 1]
    x_sq = [0, 0, 1]
    inner = add(mul(L, mul(H, H, p), p),
                scale(mul([b], mul(one_x_sq,
                                   add(mul([w], L, p),
                                       scale(x_sq, -1, p), p), p), p),
                      4, p), p)
    return mul(L, inner, p)


def shift_scale(f, c, d, p):
    # Coefficients of F(d*Z-c).
    g = [0] * 6
    for j, fj in enumerate(f):
        for i in range(j + 1):
            g[i] = (g[i] + fj * comb(j, i) * pow(d, i, p)
                    * pow(-c, j - i, p)) % p
    return g


def count_prime(p, sample_limit=8):
    if p in (2, 5):
        raise ValueError("p must be an odd prime different from 5")
    squares = {x * x % p for x in range(1, p)}
    quotient = 0
    exceptional = 0
    g0_zero_tested = 0
    samples = []
    exceptional_samples = []
    for b in range(p):
        if b == 0 or b == 1 or 2 * b % p == 1:
            continue
        for w in range(1, p):
            f = compact_f(b, w, p)
            if len(f) != 6:
                continue
            deriv = [(i * f[i]) % p for i in range(1, 6)]
            if not coprime(f, deriv, p):
                continue
            for c in range(p):
                for d in range(1, p):
                    g = shift_scale(f, c, d, p)
                    g0, g1, g2, g3, g4, g5 = g
                    if g0 == 0:
                        g0_zero_tested += 1
                        if g1 == 0 or g1 != g5:
                            continue
                        if p == 19:
                            es = range(p)
                        else:
                            es = [((6 * g1 + 8 * g4)
                                   * pow(19 * g1, -1, p)) % p]
                        for e in es:
                            if (e * e - 4) % p == 0:
                                continue
                            delta = pow(2 - e, 3, p)
                            l3 = (-g1) % p
                            l4 = (-g2) % p
                            l5 = (-g3) % p
                            l6 = (-g4) % p
                            l7 = (-g5) % p
                            c1 = (l3 - l7) % p
                            c2 = ((32 * e * e - 33 * e + 58) * l3
                                  - 20 * l5) % p
                            c3 = ((19 * e - 6) * l3 - 8 * l6) % p
                            c4 = (4 * (19 * e - 6) * l3 * l3
                                  - 32 * l4 * l3) % p
                            if c1 or c2 or c3 or c4:
                                continue
                            if not coprime([1, e, 1], g, p):
                                continue
                            tau = (-5 * delta * pow(8 * l3, -1, p)) % p
                            if tau not in squares:
                                continue
                            exceptional += 1
                            if len(exceptional_samples) < sample_limit:
                                exceptional_samples.append(
                                    (b, w, c, d, e, tau))
                        continue
                    h = (g0 + g1 - g5) % p
                    ell = (g0 - g5) % p
                    if ell == 0 or h == 5 * g0 % p or h == -5 * g0 % p:
                        continue
                    k3 = ((38 * h - 30 * g0) * ell
                          - 40 * g0 * (h - g4)) % p
                    if k3:
                        continue
                    k2 = ((128 * h * h - 330 * h * g0 + 1450 * g0 * g0)
                          * ell - 1250 * g0**3 - 150 * g0 * h * h
                          + 500 * g0 * g0 * g3) % p
                    if k2:
                        continue
                    k4 = ((95 * h - 75 * g0) * ell * ell
                          - 10 * (25 * g0 * g0 + 3 * h * h) * ell
                          + 100 * g0 * g2 * ell + (5 * g0 - h)**3) % p
                    if k4:
                        continue
                    e = 2 * h * pow(5 * g0, -1, p) % p
                    q = [1, e, 1]
                    if not coprime(q, g, p):
                        continue
                    tau = (-(5 * g0 - h)**3
                           * pow(25 * g0**3 * ell, -1, p)) % p
                    if tau not in squares:
                        continue
                    quotient += 1
                    if len(samples) < sample_limit:
                        samples.append((b, w, c, d, e, tau))
    print(f"p={p} generic_quotient={quotient} "
          f"exceptional_quotient={exceptional} "
          f"signed={2*(quotient+exceptional)} "
          f"g0_zero_tested={g0_zero_tested} samples={samples} "
          f"exceptional_samples={exceptional_samples}")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("primes", nargs="+", type=int)
    parser.add_argument("--sample-limit", type=int, default=8)
    args = parser.parse_args()
    for prime in args.primes:
        count_prime(prime, args.sample_limit)
