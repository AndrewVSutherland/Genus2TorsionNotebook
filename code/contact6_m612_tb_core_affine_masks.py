#!/usr/bin/env python3
"""Exact good-affine finite masks for the TB-half + [3,3] fiber product.

The finite enumeration builds in K=m^2 and M=L^2 and enforces every open
factor used by the rational recovery.  It is therefore an exact mask on the
good affine chart, but not a projective exclusion: a failed residue is
boundary-deferred.  When filtering rational (s,v), a prime dividing either
base denominator is retained automatically.
"""

from argparse import ArgumentParser
from fractions import Fraction
from math import gcd


BASE_PRIORITY = {
    5: {(0, 0), (0, 2), (0, 3)},
    7: {(0, v) for v in range(7)} | {(2, 3), (2, 5), (2, 6)},
    11: {(0, v) for v in range(11)}
    | {
        (1, 3), (1, 8), (3, 10), (6, 0), (6, 5), (6, 9), (6, 10),
        (7, 2), (7, 5), (7, 8), (8, 10), (9, 4), (9, 5), (9, 7),
    },
}


def half_records(p):
    """Return open TB-half records (a,r,m) indexed by s over F_p."""
    out = {}
    for s in range(1, p):
        ss = s * s % p
        b = (2 * ss - 3) % p
        records = []
        for r in range(1, p):
            for m in range(1, p):
                mm = m * m % p
                den = (s * (r * r - mm) - r) % p
                if den == 0:
                    continue
                a3_at_zero = (-3 + 4 * ss * r - 2 * mm) % p
                constant = (
                    s * (-3 * r * r + 4 * r - 3 * mm) - r * a3_at_zero
                ) % p
                a = -constant * pow(den, -1, p) % p
                if (a + b + 2) % p == 0:
                    continue
                a3 = (a - 3 + 4 * ss * r - 2 * mm) % p
                h2 = (
                    8
                    * ss
                    * (2 * ss * r * r + 2 * (a - 3) * r + 2 - mm * (2 * ss - 6))
                    - a3 * a3
                    - 32 * s * ss * r
                ) % p
                if h2 == 0:
                    records.append((a, r, m))
        out[s] = records
    return out


def affine_mask(p):
    """Enumerate the simultaneous open fiber product over F_p."""
    halves = half_records(p)
    allowed = set()
    for s, records in halves.items():
        b = (2 * s * s - 3) % p
        for v in range(1, p):
            found = False
            for a, _r, _m in records:
                c1 = (2 * a + 6) % p
                c2 = (a * a + 2 * b - 15) % p
                c3 = (2 * a * b + 22) % p
                c4 = (2 * a + b * b - 15) % p
                c5 = (2 * b + 6) % p
                for ell in range(1, p):
                    big_m = ell * ell % p
                    for u in range(p):
                        if (u * u - 4 * v * v) % p == 0:
                            continue
                        beta = (c5 * big_m + 3 * u) % p
                        delta = (
                            4 * c4 * big_m
                            + 12 * (u * u + v * v)
                            - beta * beta
                        ) % p
                        f1 = (delta * v**3 - 4 * c1 * big_m - 12 * u * v**4) % p
                        if f1:
                            continue
                        f2 = (
                            delta * delta
                            + 64 * beta * v**3
                            - 64 * c2 * big_m
                            - 192 * (u * u * v * v + v**4)
                        ) % p
                        if f2:
                            continue
                        f3 = (
                            beta * delta
                            + 16 * v**3
                            - 8 * c3 * big_m
                            - 8 * u**3
                            - 48 * u * v * v
                        ) % p
                        if f3 == 0:
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                allowed.add((s, v))
    return allowed, sum(len(records) for records in halves.values())


def rational_parameters(height):
    return sorted(
        {
            Fraction(n, d)
            for d in range(1, height + 1)
            for n in range(-height, height + 1)
            if n and gcd(abs(n), d) == 1
        }
    )


def passes(qs, qv, p, mask):
    if qs.denominator % p == 0 or qv.denominator % p == 0:
        return True
    s = qs.numerator * pow(qs.denominator, -1, p) % p
    v = qv.numerator * pow(qv.denominator, -1, p) % p
    return (s, v) in mask


def main():
    parser = ArgumentParser()
    parser.add_argument("--primes", nargs="+", type=int, default=[13, 17, 19, 23])
    parser.add_argument("--height", type=int, default=10)
    parser.add_argument("--print-masks", action="store_true")
    args = parser.parse_args()

    masks = dict(BASE_PRIORITY)
    for p in args.primes:
        mask, half_count = affine_mask(p)
        masks[p] = mask
        print(f"p {p} half_tuples {half_count} affine_allowed {len(mask)}")
        if args.print_masks:
            print("MASK", p, sorted(mask))

    values = rational_parameters(args.height)
    survivors = [(s, v) for s in values for v in values]
    print("height", args.height, "initial", len(survivors))
    for p in sorted(masks):
        survivors = [(s, v) for s, v in survivors if passes(s, v, p, masks[p])]
        print("after", p, len(survivors))
    print("SURVIVORS", survivors)


if __name__ == "__main__":
    main()
