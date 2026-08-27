#!/usr/bin/env python3
"""Necessary order/2-rank prefilter for small HLP slice directions.

This is a dependency-free companion to ``m612_hlp_direction_search.m``.
It computes #J(F_p) from #C(F_p) and #C(F_{p^2}) and the rational 2-rank
from the branch-factor degrees.  Its POSSIBLE masks (72 | #J and 2-rank
at least two) are still only upper bounds for the exact [6,12] masks.  They
are useful for staging the much rarer Magma invariant-factor calculations
and for independently checking all finite bad fibers.
"""

from __future__ import annotations

import argparse
import itertools
import math
from dataclasses import dataclass, field
from pathlib import Path


SEED = (187392, 0, -118767, 0, -118767, 0, 187392)
PRIME_STAGES = ((5, 7), (11, 13), (17, 19, 23))


def trim(f: list[int], p: int) -> list[int]:
    f = [a % p for a in f]
    while f and f[-1] == 0:
        f.pop()
    return f


def divrem(a: list[int], b: list[int], p: int):
    a, b = trim(a, p), trim(b, p)
    ib = pow(b[-1], -1, p)
    q = [0] * max(1, len(a) - len(b) + 1)
    while len(a) >= len(b):
        c = a[-1] * ib % p
        j = len(a) - len(b)
        q[j] = c
        for i, bi in enumerate(b):
            a[i + j] = (a[i + j] - c * bi) % p
        a = trim(a, p)
    return trim(q, p), a


def remainder(a: list[int], b: list[int], p: int) -> list[int]:
    return divrem(a, b, p)[1]


def monic(f: list[int], p: int) -> list[int]:
    f = trim(f, p)
    inverse = pow(f[-1], -1, p)
    return trim([inverse * c for c in f], p)


def gcd_poly(a: list[int], b: list[int], p: int) -> list[int]:
    a, b = trim(a, p), trim(b, p)
    while b:
        a, b = b, remainder(a, b, p)
    return monic(a, p)


def mulmod(a: list[int], b: list[int], modulus: list[int], p: int):
    c = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            c[i + j] = (c[i + j] + ai * bj) % p
    return remainder(c, modulus, p)


def powmod(a: list[int], exponent: int, modulus: list[int], p: int):
    result = [1]
    while exponent:
        if exponent & 1:
            result = mulmod(result, a, modulus, p)
        a = mulmod(a, a, modulus, p)
        exponent >>= 1
    return result


def sub_poly(a: list[int], b: list[int], p: int):
    return trim([
        (a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)
        for i in range(max(len(a), len(b)))
    ], p)


def factor_degrees_squarefree(f: list[int], p: int) -> list[int]:
    """Distinct-degree factorization, returning degrees with multiplicity."""
    remaining = monic(f, p)
    h = [0, 1]
    xpoly = [0, 1]
    degrees: list[int] = []
    d = 1
    while 2 * d <= len(remaining) - 1:
        h = powmod(h, p, remaining, p)
        block = gcd_poly(sub_poly(h, xpoly, p), remaining, p)
        if len(block) > 1:
            assert (len(block) - 1) % d == 0
            degrees.extend([d] * ((len(block) - 1) // d))
            remaining, rem = divrem(remaining, block, p)
            assert not rem
            remaining = monic(remaining, p)
            if len(remaining) == 1:
                break
            h = remainder(h, remaining, p)
        d += 1
    if len(remaining) > 1:
        degrees.append(len(remaining) - 1)
    return degrees


def two_rank(f: list[int], p: int) -> int:
    f = trim(f, p)
    degrees = factor_degrees_squarefree(f, p)
    s = len(degrees)
    if len(f) - 1 == 5:
        return s - 1
    return s - 1 if all(d % 2 == 0 for d in degrees) else s - 2


def is_smooth_genus2(f: list[int], p: int) -> bool:
    f = trim(f, p)
    if len(f) - 1 not in (5, 6):
        return False
    a = f
    b = trim([i * f[i] for i in range(1, len(f))], p)
    while b:
        a, b = b, remainder(a, b, p)
    return len(a) == 1


def chi(a: int, p: int) -> int:
    a %= p
    if not a:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def nonsquare(p: int) -> int:
    return next(a for a in range(2, p) if chi(a, p) == -1)


def mul2(z: tuple[int, int], w: tuple[int, int], p: int, d: int):
    """Multiply in F_p[alpha]/(alpha^2-d), with d nonsquare."""
    a, b = z
    c, e = w
    return (a * c + b * e * d) % p, (a * e + b * c) % p


def eval2(f: list[int], z: tuple[int, int], p: int, d: int):
    v = (0, 0)
    for c in reversed(f):
        v = mul2(v, z, p, d)
        v = ((v[0] + c) % p, v[1])
    return v


def jacobian_order(f: list[int], p: int) -> int:
    f = trim(f, p)
    degree = len(f) - 1
    n1 = 0
    for x in range(p):
        value = sum(c * pow(x, i, p) for i, c in enumerate(f)) % p
        n1 += 1 + chi(value, p)
    if degree == 5:
        n1 += 1
    elif chi(f[-1], p) == 1:
        n1 += 2

    # Every nonzero F_p scalar is a square in F_{p^2}, so a degree-six
    # model always has two F_{p^2}-points at infinity.
    d = nonsquare(p)
    n2 = 0
    for a in range(p):
        for b in range(p):
            u, v = eval2(f, (a, b), p, d)
            norm = (u * u - d * v * v) % p
            n2 += 1 + chi(norm, p)
    n2 += 1 if degree == 5 else 2

    a1 = p + 1 - n1
    a2 = (n2 - p * p - 1 + a1 * a1) // 2
    return 1 - a1 + a2 - p * a1 + p * p


def primitive_transverse_directions(box: int):
    for g in itertools.product(range(-box, box + 1), repeat=7):
        if not any(g):
            continue
        if math.gcd(*map(abs, g)) != 1:
            continue
        if next(a for a in g if a) < 0:
            continue
        nminus = 1298 * g[1] + 2423 * g[3] + 1298 * g[5]
        nplus = -649 * g[0] - 3072 * g[2] + 3072 * g[4] + 649 * g[6]
        if nminus and nplus:
            yield g


def normalized_mod_p(g: tuple[int, ...], p: int):
    w = tuple(a % p for a in g)
    scale = next(a for a in w if a)
    inverse = pow(scale, -1, p)
    return tuple(a * inverse % p for a in w), scale


def base_mask(n: tuple[int, ...], p: int):
    possible, bad = [], []
    for u in range(p):
        f = [(SEED[i] + u * n[i]) % p for i in range(7)]
        if not is_smooth_genus2(f, p):
            bad.append(u)
        elif jacobian_order(f, p) % 72 == 0 and two_rank(f, p) >= 2:
            possible.append(u)
    return tuple(possible), tuple(bad)


@dataclass
class Record:
    g: tuple[int, ...]
    masks: list[tuple[int, ...]] = field(default_factory=list)
    bads: list[tuple[int, ...]] = field(default_factory=list)
    primes: list[int] = field(default_factory=list)

    def key(self):
        extras = [sum(t != 0 for t in mask) for mask in self.masks]
        goods = [
            (p - 1) - sum(t != 0 for t in bad)
            for p, bad in zip(self.primes, self.bads)
        ]
        breadth = sum(n > 0 for n in extras)
        density = sum(n / good for n, good in zip(extras, goods) if good)
        return breadth, sum(extras), density, -sum(map(len, self.bads)), self.g


def add_prime(records: list[Record], p: int, emit) -> None:
    cache: dict[tuple[int, ...], tuple[tuple[int, ...], tuple[int, ...]]] = {}
    for j, record in enumerate(records, 1):
        n, scale = normalized_mod_p(record.g, p)
        if n not in cache:
            cache[n] = base_mask(n, p)
        possible, bad = cache[n]
        inverse = pow(scale, -1, p)
        record.masks.append(tuple(sorted(u * inverse % p for u in possible)))
        record.bads.append(tuple(sorted(u * inverse % p for u in bad)))
        record.primes.append(p)
        if j % 1000 == 0:
            emit(f"PRIME_PROGRESS p={p} records={j}/{len(records)} cache={len(cache)}")
    emit(f"PRIME_DONE p={p} records={len(records)} unique_projective={len(cache)}")


def report(label: str, records: list[Record], n: int, emit) -> None:
    emit(f"REPORT {label} records={len(records)} showing={min(n, len(records))}")
    for rank, r in enumerate(records[:n], 1):
        emit(f"RANK {rank} G={r.g} score={r.key()[:-1]}")
        for p, possible, bad in zip(r.primes, r.masks, r.bads):
            emit(f" MASK p={p} possible72_r2={possible} bad={bad} infinity=SEPARATE")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--box", type=int, default=2)
    ap.add_argument("--stage-top", type=int, default=5000)
    ap.add_argument("--final-top", type=int, default=250)
    ap.add_argument("--report-top", type=int, default=20)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    lines: list[str] = []

    def emit(line: str) -> None:
        print(line, flush=True)
        lines.append(line)

    records = [Record(g) for g in primitive_transverse_directions(args.box)]
    emit(f"M612_HLP_DIRECTION_ORDER_PREFILTER box={args.box} directions={len(records)}")
    for p in PRIME_STAGES[0]:
        add_prime(records, p, emit)
    records.sort(key=Record.key, reverse=True)
    report("STAGE1", records, args.report_top, emit)
    records = records[: args.stage_top]

    for p in PRIME_STAGES[1]:
        add_prime(records, p, emit)
    records.sort(key=Record.key, reverse=True)
    report("STAGE2", records, args.report_top, emit)
    records = records[: args.final_top]

    for p in PRIME_STAGES[2]:
        add_prime(records, p, emit)
    records.sort(key=Record.key, reverse=True)
    report("FINAL", records, args.report_top, emit)

    baseline = [Record((1, 1, 0, 0, 0, 0, 0))]
    for stage in PRIME_STAGES:
        for p in stage:
            add_prime(baseline, p, lambda _: None)
    report("BASELINE_1_PLUS_X", baseline, 1, emit)
    emit("M612_HLP_DIRECTION_ORDER_PREFILTER_DONE")

    if args.output:
        args.output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
