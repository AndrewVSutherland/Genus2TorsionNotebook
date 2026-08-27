#!/usr/bin/env python3
"""Projective CRT height sieve for the transverse HLP marked slice."""

from __future__ import annotations

import argparse
from math import gcd
from pathlib import Path


def parse_masks(paths: list[Path]) -> dict[int, tuple[int, ...]]:
    masks = {}
    for path in paths:
        for line in path.read_text().splitlines():
            if not line or line.startswith("#"):
                continue
            left, rest = line.split(":", 1)
            p = int(left.strip())
            field = rest.split(";", 1)[0].strip()
            masks[p] = tuple(int(x) for x in field.split(",") if x)
    return masks


def invmod(a: int, p: int) -> int:
    return pow(a % p, -1, p)


def crt_extend(residue: int, modulus: int, target: int, p: int) -> int:
    step = ((target - residue) % p) * invmod(modulus, p) % p
    return residue + modulus * step


def ceil_div(a: int, b: int) -> int:
    return -((-a) // b)


def base_residues(d: int, base: list[int], masks: dict[int, tuple[int, ...]]):
    states = [(0, 1)]
    for p in base:
        next_states = []
        dm = d % p
        targets = range(1, p) if dm == 0 else \
            ((dm * a) % p for a in masks[p])
        targets = tuple(targets)
        for residue, modulus in states:
            for target in targets:
                next_states.append((crt_extend(residue, modulus, target, p),
                                    modulus * p))
        states = next_states
    return states


def passes(n: int, d: int, tests: list[int], masks: dict[int, tuple[int, ...]]):
    for p in tests:
        dm = d % p
        if dm == 0:
            continue  # projective t=infinity is conservatively allowed
        t = (n % p) * invmod(dm, p) % p
        if t not in masks[p]:
            return False, p
    return True, 0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--height", type=int, default=100000)
    ap.add_argument("--masks", type=Path, nargs="+", required=True)
    ap.add_argument("--base", default="11,13,31,43")
    ap.add_argument("--out", type=Path)
    ap.add_argument("--progress", type=int, default=10000)
    args = ap.parse_args()

    masks = parse_masks(args.masks)
    base = [int(x) for x in args.base.split(",") if x]
    if any(p not in masks for p in base):
        raise ValueError("base prime missing from masks")
    tests = sorted((p for p in masks if p not in base),
                   key=lambda p: (len(masks[p]) + 1) / (p + 1))
    H = args.height
    survivors = []
    leaves = bounded = primitive = 0
    kill = {p: 0 for p in tests}

    print("HLP transverse projective CRT sieve")
    print("height", H, "base", base,
          "base_modulus", __import__("math").prod(base), "tests", tests)

    for d in range(1, H + 1):
        for residue, modulus in base_residues(d, base, masks):
            leaves += 1
            k0 = ceil_div(-H - residue, modulus)
            k1 = (H - residue) // modulus
            for k in range(k0, k1 + 1):
                n = residue + k * modulus
                bounded += 1
                if gcd(abs(n), d) != 1:
                    continue
                primitive += 1
                ok, badp = passes(n, d, tests, masks)
                if not ok:
                    kill[badp] += 1
                    continue
                survivors.append((n, d))
                print(f"SURVIVOR n={n} d={d} t={n}/{d}")
        if args.progress and d % args.progress == 0:
            print("progress d", d, "leaves", leaves, "primitive", primitive,
                  "survivors", len(survivors))

    print("DONE height", H, "leaves", leaves, "bounded", bounded,
          "primitive", primitive, "survivors", len(survivors))
    print("kill_counts", [(p, kill[p]) for p in tests])
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        lines = ["# n d t"] + [f"{n} {d} {n}/{d}" for n, d in survivors]
        args.out.write_text("\n".join(lines) + "\n")
        print("wrote", args.out)


if __name__ == "__main__":
    main()
