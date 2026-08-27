#!/usr/bin/env python3
"""Count effective t=eps*s mod 7 charts in the exact height-80 box."""

import argparse
from collections import Counter
from fractions import Fraction
from math import gcd


def rational_values(height):
    return sorted({
        Fraction(numerator, denominator)
        for denominator in range(1, height + 1)
        for numerator in range(-height, height + 1)
        if gcd(numerator, denominator) == 1
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=80)
    args = parser.parse_args()

    counts = Counter()
    excluded = Counter()
    values = rational_values(args.height)
    for s in values:
        for epsilon in (-1, 1):
            t = epsilon * s
            # These are exactly the eight nonsmooth/undefined cases in the
            # recorded height-80 Magma run.
            if t in (Fraction(-1), Fraction(1), Fraction(0), Fraction(-1, 2)):
                excluded[str(t)] += 1
                continue
            if t.denominator % 7 == 0:
                counts["denominator_divisible_7"] += 1
                continue
            residue = t.numerator * pow(t.denominator, -1, 7) % 7
            counts[f"t_mod7_{residue}"] += 1

    print(f"height={args.height}")
    print(f"rational_s_count={len(values)}")
    print(f"smooth_marked9_count={sum(counts.values())}")
    print(f"excluded_count={sum(excluded.values())}")
    print(f"excluded={sorted(excluded.items())}")
    print(f"chart_counts={sorted(counts.items())}")
    print("chart_meanings:")
    print("  t_mod7_2 = only smooth displayed mod-7 chart; first-killed at p=7")
    print("  t_mod7_4 = live ordinary-node q-cover chart")
    print("  t_mod7_6 = t=-1 pole chart")
    print("  denominator_divisible_7 = no direct reduction of the displayed t-chart")


if __name__ == "__main__":
    main()
