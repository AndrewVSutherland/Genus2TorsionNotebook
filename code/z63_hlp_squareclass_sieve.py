#!/usr/bin/env python3
"""Low-memory necessary sieve for HLP 7-by-9 gluings.

For reduced rationals t,u of projective height at most H, group the values

  Delta_7(t) = t(t-1)(t^3-8t^2+5t+1),
  Delta_9(u) = u(u-1)(u^2-u+1)(u^3-6u^2+3u+1)

by their exact rational squareclass.  Equal squareclass is necessary, but not
sufficient, for the two 2-torsion cubic fields to be isomorphic.  The output
is deliberately compact; candidate pairs are passed to Magma for the cubic-
field isomorphism test.
"""

from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction
from math import gcd

from sympy import factorint


def rationals(height):
    for den in range(1, height + 1):
        for num in range(-height, height + 1):
            if gcd(abs(num), den) == 1:
                yield Fraction(num, den)


def delta7(t):
    return t * (t - 1) * (t**3 - 8 * t**2 + 5 * t + 1)


def delta9(u):
    return u * (u - 1) * (u**2 - u + 1) * (
        u**3 - 6 * u**2 + 3 * u + 1
    )


def squareclass(q):
    """Canonical squarefree integer representing q in Q*/Q*2."""
    n = q.numerator * q.denominator
    sign = -1 if n < 0 else 1
    out = sign
    for p, exponent in factorint(abs(n)).items():
        if exponent & 1:
            out *= int(p)
    return out


def main():
    parser = ArgumentParser()
    parser.add_argument("--height", type=int, default=20)
    parser.add_argument("--show", type=int, default=200,
                        help="maximum candidate pairs to print")
    args = parser.parse_args()

    values = list(rationals(args.height))
    by7 = defaultdict(list)
    by9 = defaultdict(list)
    for t in values:
        d = delta7(t)
        if d:
            by7[squareclass(d)].append(t)
    for u in values:
        d = delta9(u)
        if d:
            by9[squareclass(d)].append(u)

    common = sorted(set(by7).intersection(by9))
    pairs = [(t, u, sc) for sc in common for t in by7[sc] for u in by9[sc]]
    known = (Fraction(-16, 3), Fraction(4, 1))

    print(f"Z63_HLP_SQUARECLASS_SIEVE height={args.height}")
    print(f"parameters_each={len(values)}")
    print(f"squareclasses_7={len(by7)} squareclasses_9={len(by9)}")
    print(f"common_squareclasses={len(common)} candidate_pairs={len(pairs)}")
    print(f"known_pair_present={any((t,u)==known for t,u,_ in pairs)}")
    for t, u, sc in pairs[: args.show]:
        print(f"PAIR t={t} u={u} squareclass={sc}")
    if len(pairs) > args.show:
        print(f"TRUNCATED omitted={len(pairs)-args.show}")
    print("Z63_HLP_SQUARECLASS_SIEVE_DONE")


if __name__ == "__main__":
    main()
