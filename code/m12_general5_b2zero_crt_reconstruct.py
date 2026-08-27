#!/usr/bin/env python3
"""CRT/rational reconstruction on the generic linear-B M(12)+5 quotient.

The quotient curve has very few fully open points over the small finite
fields.  This script enumerates all points whose sign-cover coordinate tau
is a nonzero square, combines their coordinates by CRT, and applies the
canonical balanced rational reconstruction.  Every reconstructed tuple is
then substituted in the exact characteristic-zero equations; no match is
reported from CRT compatibility alone.

By default we test all subsets containing at least five of
7,11,13,17,19,23,29.  Omitting one or two primes permits a rational point to
have bad/boundary reduction there, while retaining a useful reconstruction
modulus.  Memory consumption is O(sum #points), not O(product p).
"""

from argparse import ArgumentParser
from fractions import Fraction
from itertools import combinations, product
from math import comb, gcd, isqrt, prod

from m12_general5_b2zero_rootquotient_local import (
    compact_f,
    coprime,
    shift_scale,
)


def signed_quotient_points(p):
    """Return fully open (b,w,c,d,e,tau) quotient points over F_p."""
    if p in (2, 5):
        raise ValueError("p must be an odd prime other than 5")
    squares = {x * x % p for x in range(1, p)}
    ans = []
    for b in range(p):
        if b in (0, 1) or 2 * b % p == 1:
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
                    if g0:
                        h = (g0 + g1 - g5) % p
                        ell = (g0 - g5) % p
                        if (ell == 0 or h == 5 * g0 % p
                                or h == -5 * g0 % p):
                            continue
                        k3 = ((38 * h - 30 * g0) * ell
                              - 40 * g0 * (h - g4)) % p
                        if k3:
                            continue
                        k2 = ((128 * h * h - 330 * h * g0
                               + 1450 * g0 * g0) * ell
                              - 1250 * g0**3 - 150 * g0 * h * h
                              + 500 * g0 * g0 * g3) % p
                        if k2:
                            continue
                        k4 = ((95 * h - 75 * g0) * ell * ell
                              - 10 * (25 * g0 * g0 + 3 * h * h) * ell
                              + 100 * g0 * g2 * ell
                              + (5 * g0 - h)**3) % p
                        if k4:
                            continue
                        e = 2 * h * pow(5 * g0, -1, p) % p
                        if not coprime([1, e, 1], g, p):
                            continue
                        tau = (-(5 * g0 - h)**3
                               * pow(25 * g0**3 * ell, -1, p)) % p
                        if tau in squares:
                            ans.append((b, w, c, d, e, tau))
                        continue

                    # Exceptional g0=0 chart.  Work directly with the four
                    # C-equations, since 19 is exceptional in the displayed
                    # rational e formula.
                    if g1 == 0 or g1 != g5:
                        continue
                    es = range(p) if p == 19 else [
                        (6 * g1 + 8 * g4) * pow(19 * g1, -1, p) % p
                    ]
                    for e in es:
                        if (e * e - 4) % p == 0:
                            continue
                        inv2 = pow(2, -1, p)
                        inv8 = pow(8, -1, p)
                        l3 = -g1 % p
                        l4 = -g2 % p
                        l5 = -g3 % p
                        l6 = -g4 % p
                        l7 = -g5 % p
                        delta = pow(2 - e, 3, p)
                        cs = (
                            l3 - l7,
                            (32 * e * e - 33 * e + 58) * l3 - 20 * l5,
                            (19 * e - 6) * l3 - 8 * l6,
                            4 * (19 * e - 6) * l3 * l3
                            - 32 * ((5 * inv2 + 15 * e * e * inv8)
                                    * g0 - g2) * l3
                            + 5 * delta * g0 * g0,
                        )
                        if any(x % p for x in cs):
                            continue
                        if not coprime([1, e, 1], g, p):
                            continue
                        tau = -5 * delta * pow(8 * l3, -1, p) % p
                        if tau in squares:
                            ans.append((b, w, c, d, e, tau))
    return ans


def crt_pair(a, m, b, n):
    """CRT for coprime positive moduli, returning a residue modulo mn."""
    return (a + ((b - a) * pow(m, -1, n) % n) * m) % (m * n)


def crt_coordinate(residues, primes):
    a, modulus = 0, 1
    for b, p in zip(residues, primes):
        a = crt_pair(a, modulus, b, p)
        modulus *= p
    return a, modulus


def balanced_reconstruction(a, modulus):
    """Canonical reconstruction with |num|,den <= floor(sqrt(m/2))."""
    bound = isqrt((modulus - 1) // 2)
    if a == 0:
        return Fraction(0), bound
    r0, r1 = modulus, a % modulus
    t0, t1 = 0, 1
    while r1 > bound:
        q = r0 // r1
        r0, r1 = r1, r0 - q * r1
        t0, t1 = t1, t0 - q * t1
    num, den = r1, t1
    if den < 0:
        num, den = -num, -den
    if (den == 0 or den > bound or abs(num) > bound
            or gcd(abs(num), den) != 1
            or (num - a * den) % modulus):
        return None, bound
    return Fraction(num, den), bound


def qadd(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0)
            + (b[i] if i < len(b) else 0) for i in range(n)]


def qscale(a, scalar):
    return [scalar * x for x in a]


def qmul(a, b):
    ans = [Fraction(0)] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            ans[i + j] += x * y
    return ans


def compact_f_q(b, w):
    L = [b, 2 * b - 1]
    H = [w, 1 + w * b]
    inner = qadd(
        qmul(L, qmul(H, H)),
        qscale(qmul([b], qmul(
            [1, 2, 1],
            qadd(qmul([w], L), [0, 0, -1]))), 4))
    return qmul(L, inner)


def shift_scale_q(f, c, d):
    g = [Fraction(0)] * 6
    for j, fj in enumerate(f):
        for i in range(j + 1):
            g[i] += fj * comb(j, i) * d**i * (-c)**(j - i)
    return g


def is_square_q(x):
    if x <= 0:
        return False
    return (isqrt(x.numerator)**2 == x.numerator
            and isqrt(x.denominator)**2 == x.denominator)


def exact_quotient_check(point):
    b, w, c, d, e, tau_reconstructed = point
    if (b == 0 or w == 0 or b == 1 or 2 * b == 1 or d == 0
            or e in (2, -2)):
        return False
    f = compact_f_q(b, w)
    g0, g1, g2, g3, g4, g5 = shift_scale_q(f, c, d)
    ell3 = Fraction(5, 2) * e * g0 - g1
    ell4 = (Fraction(5, 2) + Fraction(15, 8) * e * e) * g0 - g2
    ell5 = (Fraction(5, 2) + Fraction(15, 8) * e * e) * g0 - g3
    ell6 = Fraction(5, 2) * e * g0 - g4
    ell7 = g0 - g5
    delta = (2 - e)**3
    cs = (
        ell3 - ell7,
        (32 * e * e - 33 * e + 58) * ell3 - 20 * ell5,
        (19 * e - 6) * ell3 - 8 * ell6,
        4 * (19 * e - 6) * ell3 * ell3 - 32 * ell4 * ell3
        + 5 * delta * g0 * g0,
    )
    if any(cs) or ell3 == 0:
        return False
    tau = -5 * delta / (8 * ell3)
    return tau == tau_reconstructed and is_square_q(tau)


def main():
    parser = ArgumentParser()
    parser.add_argument("--primes", nargs="+", type=int,
                        default=[7, 11, 13, 17, 19, 23, 29])
    parser.add_argument("--min-subset", type=int, default=5)
    args = parser.parse_args()
    primes = tuple(args.primes)
    local = {}
    for p in primes:
        local[p] = signed_quotient_points(p)
        print(f"LOCAL p={p} quotient_square_tau={len(local[p])} "
              f"points={local[p]}")

    tested = reconstructed = 0
    exact = set()
    for size in range(len(primes), args.min_subset - 1, -1):
        for subset in combinations(primes, size):
            modulus = prod(subset)
            bound = isqrt((modulus - 1) // 2)
            combos = prod(len(local[p]) for p in subset)
            subset_reconstructed = 0
            for choices in product(*(local[p] for p in subset)):
                tested += 1
                values = []
                for coordinate in range(6):
                    residue, m = crt_coordinate(
                        [point[coordinate] for point in choices], subset)
                    value, _ = balanced_reconstruction(residue, m)
                    if value is None:
                        break
                    values.append(value)
                if len(values) != 6:
                    continue
                reconstructed += 1
                subset_reconstructed += 1
                candidate = tuple(values)
                if exact_quotient_check(candidate):
                    exact.add(candidate)
                    print(f"EXACT_CANDIDATE subset={subset} point={candidate}")
            print(f"SUBSET primes={subset} modulus={modulus} bound={bound} "
                  f"combinations={combos} reconstructed={subset_reconstructed}")
    print(f"DONE tested={tested} reconstructed={reconstructed} "
          f"exact_candidates={len(exact)}")


if __name__ == "__main__":
    main()
