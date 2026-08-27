#!/usr/bin/env python3
"""Meet-in-the-middle search on the exact orbit-01 halving cover.

After scaling r1=1, write

    r_j = (1-t_j^2)/(1+t_j^2),
    z_j = 2*t_j/(1+t_j^2).

Then r1^2-r_j^2=z_j^2 automatically, and the two Clebsch--Klein
equations are exactly

    sum_j 1/(1+t_j^2) = 3/2,
    product_j (1-t_j^2) = 16.

The equations split into pair invariants, so a hash join replaces an
O(N^4) search by O(N^2).  Height means max(|numerator|,denominator)
for each nonnegative t; signs of t do not change r.
"""

from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from math import gcd
from time import monotonic


def values(height):
    out = []
    for d in range(1, height + 1):
        for n in range(1, height + 1):
            if gcd(n, d) != 1 or n == d:
                continue
            t = Fraction(n, d)
            a = 1 / (1 + t * t)
            b = 1 - t * t
            r = 2 * a - 1
            out.append((t, a, b, r))
    out.sort(key=lambda row: row[0])
    return out


def main():
    ap = ArgumentParser()
    ap.add_argument("height", nargs="?", type=int, default=30)
    ap.add_argument("--max-hits", type=int, default=20)
    args = ap.parse_args()
    vals = values(args.height)
    start = monotonic()
    table = defaultdict(list)
    pairs = 0
    joins = 0
    hits = set()

    # Store unordered pairs.  t and 1/t give opposite r and hence the
    # same branch square, so exclude them within a prospective curve.
    for i in range(len(vals)):
        ti, ai, bi, ri = vals[i]
        for j in range(i + 1, len(vals)):
            tj, aj, bj, rj = vals[j]
            if ri * ri == rj * rj:
                continue
            pairs += 1
            key = (ai + aj, bi * bj)
            want_s = Fraction(3, 2) - key[0]
            if key[1] != 0:
                want = (want_s, Fraction(16, 1) / key[1])
                for k, ell in table.get(want, ()):
                    joins += 1
                    inds = (k, ell, i, j)
                    rs = [vals[q][3] for q in inds]
                    if len({r * r for r in rs}) != 4 or any(r == 0 for r in rs):
                        continue
                    ts = tuple(sorted(vals[q][0] for q in inds))
                    # Exact independent verification of both equations.
                    if sum((1 / (1 + t * t) for t in ts), Fraction()) != Fraction(3, 2):
                        raise AssertionError("sum identity failed")
                    prod = Fraction(1)
                    for t in ts:
                        prod *= 1 - t * t
                    if prod != 16:
                        raise AssertionError("product identity failed")
                    hits.add(ts)
                    if len(hits) >= args.max_hits:
                        break
            table[key].append((i, j))
        if len(hits) >= args.max_hits:
            break

    print("ELKIES22210_ORBIT01_UNIT_CIRCLE_SEARCH")
    print("height", args.height)
    print("values", len(vals))
    print("pairs", pairs)
    print("pair_keys", len(table))
    print("joins", joins)
    print("hits", len(hits))
    for ts in sorted(hits):
        rs = tuple((1 - t * t) / (1 + t * t) for t in ts)
        print("HIT t", ts, "r", (Fraction(1),) + rs)
    print("elapsed_seconds", f"{monotonic()-start:.3f}")
    print("DONE")


if __name__ == "__main__":
    main()
