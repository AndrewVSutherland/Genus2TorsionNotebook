#!/usr/bin/env python3
"""Projective local lift check for the M(2,2,4,8) q^2=1 diagonal.

On the q^2=1 diagonal of the D-square + F0-square surface, the independent
full-cover square classes reduce to F1 and F4.  After homogenizing
d=D/T, n=N/T and clearing square denominators, necessary local conditions are

    H1 = -8*(D^2-T^2)*(D-N+T)*(D+N-T) in Q_p^{*2}
    H4 =  8*G_h(D,N,T)                         in Q_p^{*2}.

The earlier affine-good check only looked at T=1 with all displayed
denominators units.  This script searches all three projective unit charts and
lifts residues modulo p^k.  Empty survivors at some depth are a genuine
projective local obstruction for this diagonal.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from typing import Iterable


@dataclass(frozen=True)
class Chart:
    name: str
    fixed: str

    def point(self, a: int, b: int) -> tuple[int, int, int]:
        if self.fixed == "T":
            return a, b, 1
        if self.fixed == "D":
            return 1, a, b
        if self.fixed == "N":
            return a, 1, b
        raise ValueError(self.fixed)


CHARTS = [Chart("T=1", "T"), Chart("D=1", "D"), Chart("N=1", "N")]


def h1_value(D: int, N: int, T: int) -> int:
    return -8 * (D * D - T * T) * (D - N + T) * (D + N - T)


def gh_value(D: int, N: int, T: int) -> int:
    D2 = D * D
    N2 = N * N
    T2 = T * T
    return (
        D2 * D2
        + D2 * N2
        + 2 * D2 * N * T
        - 2 * D2 * T2
        + N2 * T2
        - 2 * N * T2 * T
        + T2 * T2
    )


def h4_value(D: int, N: int, T: int) -> int:
    return 8 * gh_value(D, N, T)


def square_residues(modulus: int) -> set[int]:
    return {(x * x) % modulus for x in range(modulus)}


def keep_pair(chart: Chart, a: int, b: int, modulus: int, squares: set[int]) -> bool:
    D, N, T = chart.point(a, b)
    return (h1_value(D, N, T) % modulus in squares) and (
        h4_value(D, N, T) % modulus in squares
    )


def resolved_nonzero(chart: Chart, a: int, b: int, modulus: int) -> bool:
    D, N, T = chart.point(a, b)
    return h1_value(D, N, T) % modulus != 0 and h4_value(D, N, T) % modulus != 0


def unit_nonzero(chart: Chart, a: int, b: int, p: int) -> bool:
    D, N, T = chart.point(a, b)
    return h1_value(D, N, T) % p != 0 and h4_value(D, N, T) % p != 0


def lift_chart(chart: Chart, p: int, depth: int) -> list[list[tuple[int, int]]]:
    survivors_by_depth: list[list[tuple[int, int]]] = []
    modulus = p
    squares = square_residues(modulus)
    survivors = [
        (a, b)
        for a in range(modulus)
        for b in range(modulus)
        if keep_pair(chart, a, b, modulus, squares)
    ]
    survivors_by_depth.append(survivors)

    for k in range(1, depth):
        old_modulus = modulus
        modulus *= p
        squares = square_residues(modulus)
        lifted: list[tuple[int, int]] = []
        for a, b in survivors:
            for da in range(p):
                aa = a + da * old_modulus
                for db in range(p):
                    bb = b + db * old_modulus
                    if keep_pair(chart, aa, bb, modulus, squares):
                        lifted.append((aa, bb))
        survivors = lifted
        survivors_by_depth.append(survivors)
        if not survivors:
            break

    return survivors_by_depth


def format_sample(chart: Chart, pair: tuple[int, int], modulus: int) -> str:
    D, N, T = chart.point(*pair)
    return (
        f"(D,N,T)=({D % modulus},{N % modulus},{T % modulus}); "
        f"H1={h1_value(D, N, T) % modulus}; H4={h4_value(D, N, T) % modulus}"
    )


def parse_primes(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part.strip()]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", default="3,5,13")
    parser.add_argument("--depth", type=int, default=6)
    parser.add_argument("--samples", type=int, default=3)
    args = parser.parse_args()

    print("AGENT M2248 q^2=1 projective local lift check")
    print("conditions H1=-8*(D^2-T^2)*(D-N+T)*(D+N-T), H4=8*G_h")
    print(f"primes {parse_primes(args.primes)} depth {args.depth}")

    for p in parse_primes(args.primes):
        print(f"PRIME_START {p}")
        active_any = False
        obstruction_depth: int | None = None
        chart_data: list[tuple[Chart, list[list[tuple[int, int]]]]] = []

        for chart in CHARTS:
            lifts = lift_chart(chart, p, args.depth)
            chart_data.append((chart, lifts))
            if lifts and lifts[-1]:
                active_any = True

        max_seen_depth = max(len(lifts) for _, lifts in chart_data)
        for k in range(1, max_seen_depth + 1):
            modulus = p**k
            counts = []
            total = 0
            resolved_total = 0
            unit_total = 0
            for chart, lifts in chart_data:
                current = lifts[k - 1] if k <= len(lifts) else []
                count = len(current)
                resolved = sum(1 for a, b in current if resolved_nonzero(chart, a, b, modulus))
                unit_count = sum(1 for a, b in current if unit_nonzero(chart, a, b, p))
                counts.append(f"{chart.name}:{count}/{resolved}/{unit_count}")
                total += count
                resolved_total += resolved
                unit_total += unit_count
            print(
                f"DEPTH p={p} k={k} mod={modulus} "
                f"total={total} resolved_nonzero={resolved_total} unit={unit_total} "
                f"{' '.join(counts)}"
            )
            if total == 0 and obstruction_depth is None:
                obstruction_depth = k

        if obstruction_depth is not None:
            print(f"PROJECTIVE_OBSTRUCTION p={p} depth={obstruction_depth}")
        else:
            print(f"NO_OBSTRUCTION_TO_DEPTH p={p} depth={args.depth}")
            final_modulus = p ** args.depth
            for chart, lifts in chart_data:
                if len(lifts) < args.depth:
                    continue
                for pair in lifts[-1][: args.samples]:
                    print(f"SURVIVOR p={p} chart={chart.name} {format_sample(chart, pair, final_modulus)}")

        if not active_any and obstruction_depth is None:
            print(f"PROJECTIVE_OBSTRUCTION p={p} depth<={args.depth}")
        print(f"PRIME_DONE {p}")

    print("DONE")


if __name__ == "__main__":
    main()
