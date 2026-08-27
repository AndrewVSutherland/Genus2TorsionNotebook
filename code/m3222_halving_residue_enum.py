#!/usr/bin/env python3
"""Fast residue enumerator for M_1(8,2,2) order-8 halving.

This is a necessary-condition sieve.  For each good open prime p, the
distinguished order-8 class must reduce to a class divisible by 2 in J(F_p).
If a rational parameter has boundary reduction at p, it is allowed through.
"""

from __future__ import annotations

import argparse
import math
from fractions import Fraction
from pathlib import Path


def read_allowed(path: Path) -> dict[int, set[int]]:
    allowed: dict[int, set[int]] = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split()
        p = int(parts[0])
        allowed[p] = {int(x) for x in parts[1:]}
    return allowed


def rational_parameters(height: int) -> list[Fraction]:
    vals: list[Fraction] = []
    seen: set[Fraction] = set()
    for den in range(1, height + 1):
        for num in range(-height, height + 1):
            if math.gcd(num, den) != 1:
                continue
            q = Fraction(num, den)
            if q not in seen:
                seen.add(q)
                vals.append(q)
    return vals


def residue(q: Fraction, p: int) -> int | None:
    den = q.denominator % p
    if den == 0:
        return None
    return (q.numerator % p) * pow(den, -1, p) % p


def open_residue(u: int, v: int, p: int) -> bool:
    s = (u + v) % p
    pp = (u * v) % p
    factors = [
        pp,
        (s + 1) % p,
        (s + 2) % p,
        (pp - s + 1) % p,
        (s - pp + 1) % p,
        (s * s - 4 * pp) % p,
        (2 * s * s + 3 * s + pp + 1) % p,
        (s * s * s - s * s * pp + s * s - 4 * s * pp - 4 * pp) % p,
    ]
    return all(x != 0 for x in factors)




def rational_open(u: Fraction, v: Fraction) -> bool:
    s = u + v
    pp = u * v
    factors = [
        pp,
        s + 1,
        s + 2,
        pp - s + 1,
        s - pp + 1,
        s * s - 4 * pp,
        2 * s * s + 3 * s + pp + 1,
        s * s * s - s * s * pp + s * s - 4 * s * pp - 4 * pp,
    ]
    return all(x != 0 for x in factors)

def passes(q1: Fraction, q2: Fraction, allowed: dict[int, set[int]]) -> tuple[bool, int]:
    for p, allowed_keys in allowed.items():
        u = residue(q1, p)
        v = residue(q2, p)
        if u is None or v is None:
            continue
        if open_residue(u, v, p):
            key = u + p * v
            if key not in allowed_keys:
                return False, p
    return True, 0


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=50)
    parser.add_argument("--allowed", type=Path, default=Path("data/m3222_halving_allowed_residues_p43.txt"))
    parser.add_argument("--out", type=Path, default=Path("data/m3222_halving_candidates_h50_p43.txt"))
    parser.add_argument("--progress", type=int, default=1_000_000)
    parser.add_argument("--require-rational-open", action="store_true")
    args = parser.parse_args()

    allowed = read_allowed(args.allowed)
    params = rational_parameters(args.height)
    args.out.parent.mkdir(parents=True, exist_ok=True)

    checked = 0
    survivors = 0
    rejected: dict[int, int] = {p: 0 for p in allowed}
    with args.out.open("w") as out:
        out.write(f"# height={args.height} primes={sorted(allowed)} params={len(params)}\n")
        out.write("# u v\n")
        for u in params:
            for v in params:
                checked += 1
                if args.require_rational_open and not rational_open(u, v):
                    continue
                ok, badp = passes(u, v, allowed)
                if not ok:
                    rejected[badp] = rejected.get(badp, 0) + 1
                    continue
                survivors += 1
                out.write(f"{u} {v}\n")
                if args.progress and checked % args.progress == 0:
                    print("progress", checked, "survivors", survivors, flush=True)

    print("DONE")
    print("height", args.height)
    print("params", len(params))
    print("checked", checked)
    print("survivors", survivors)
    print("require_rational_open", args.require_rational_open)
    for p in sorted(rejected):
        print("rejected_at", p, rejected[p])


if __name__ == "__main__":
    main()
