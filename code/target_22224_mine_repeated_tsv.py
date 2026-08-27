#!/usr/bin/env python3
"""Mine streamed A(2,2,2,8) TSVs for repeated projective triples.

If two primitive signed quadruples have the same three coordinates up to a
common rational scale and the remaining coordinates differ by a rational
square, they lie on a common one-parameter full-cover fibre

    (a,b,c,d_0 T^2).

The input may be a partially written TSV: malformed/incomplete final rows are
silently ignored, which makes this useful while a large box run is active.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path


def square_q(q: Fraction) -> bool:
    if q <= 0:
        return False
    return isqrt(q.numerator) ** 2 == q.numerator and isqrt(q.denominator) ** 2 == q.denominator


def sqrt_q(q: Fraction) -> Fraction:
    return Fraction(isqrt(q.numerator), isqrt(q.denominator))


def read_tuples(paths: list[Path]) -> set[tuple[int, int, int, int]]:
    out: set[tuple[int, int, int, int]] = set()
    for path in paths:
        if not path.exists():
            continue
        with path.open(errors="replace") as handle:
            header = handle.readline().rstrip("\n").split("\t")
            try:
                cols = [header.index(f"primitive_{x}") for x in "abcd"]
            except ValueError:
                continue
            for line in handle:
                z = line.rstrip("\n").split("\t")
                if len(z) <= max(cols):
                    continue
                try:
                    v = tuple(int(z[i]) for i in cols)
                except ValueError:
                    continue
                if len(set(v)) == 4 and all(v):
                    out.add(v)  # type: ignore[arg-type]
    return out


def charts(v: tuple[int, int, int, int]):
    for omit in range(4):
        fixed = [v[i] for i in range(4) if i != omit]
        variable = v[omit]
        h = 0
        for z in fixed:
            h = gcd(h, abs(z))
        fixed = [z // h for z in fixed]
        q = Fraction(variable, h)
        fixed.sort()
        if fixed[0] < 0:
            fixed = [-z for z in reversed(fixed)]
            q = -q
        yield tuple(fixed), q, v


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("paths", nargs="+", type=Path)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    tuples = read_tuples(args.paths)
    groups: dict[tuple[int, int, int], dict[Fraction, set[tuple[int, int, int, int]]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for v in tuples:
        for fixed, q, source in charts(v):
            groups[fixed][q].add(source)

    rows = []
    for fixed, qmap in groups.items():
        qs = sorted(qmap)
        for i, q0 in enumerate(qs):
            for q1 in qs[i + 1 :]:
                ratio = q1 / q0
                if not square_q(ratio):
                    continue
                rows.append((fixed, q0, q1, sqrt_q(ratio), qmap[q0], qmap[q1]))
    rows.sort(key=lambda r: (max(abs(z) for z in r[0]), r[0], r[1], r[2]))

    lines = [f"tuples\t{len(tuples)}", f"square_repeated_pairs\t{len(rows)}"]
    for fixed, q0, q1, t1, src0, src1 in rows:
        lines.append(
            "fixed=%s\td0=%s\td1=%s\tT1=%s\tsource0=%s\tsource1=%s"
            % (fixed, q0, q1, t1, sorted(src0), sorted(src1))
        )
    text = "\n".join(lines) + "\n"
    if args.output:
        args.output.write_text(text)
    print(text, end="")


if __name__ == "__main__":
    main()
