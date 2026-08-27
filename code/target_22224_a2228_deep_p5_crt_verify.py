#!/usr/bin/env python3
"""Independent exact verifier for deep-p^5 CRT reconstructions."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


def canonical(v: tuple[int, ...], p: int, modulus: int) -> tuple[int, ...]:
    profiles = []
    for x in v:
        if x % p:
            q = pow(x, -1, modulus)
            z = [y * q % modulus for y in v]
            z = [min(y, (-y) % modulus) for y in z]
            profiles.append(tuple(sorted(z)))
    if not profiles:
        raise ValueError("no projective unit coordinate")
    return min(profiles)


def deep_keys(path: Path, p: int, modulus: int) -> set[tuple[int, ...]]:
    with path.open() as f:
        return {
            canonical(tuple(int(r[k]) for k in "abcd"), p, modulus)
            for r in csv.DictReader(f, delimiter="\t")
        }


def primitive(v: tuple[int, ...]) -> tuple[int, ...]:
    g = math.gcd(*v)
    if g:
        v = tuple(z // abs(g) for z in v)
    for z in v:
        if z:
            return tuple(-x for x in v) if z < 0 else v
    return v


def square(n: int) -> bool:
    if n < 0:
        return False
    r = math.isqrt(n)
    return r * r == n


def radicands(v: tuple[int, int, int, int]) -> tuple[int, ...]:
    a, b, c, d = v
    return (
        a * b * c * d,
        a * (a + b) * (a + c) * (a + d),
        b * (b + a) * (b + c) * (b + d),
        c * (c + a) * (c + b) * (c + d),
    )


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", default="results/target_22224_a2228_deep_p5_crt_candidates.tsv")
    ap.add_argument("--output", default="results/target_22224_a2228_deep_p5_crt_verified.tsv")
    ap.add_argument(
        "--p11",
        default="results/target_22224_direct_contact_deep13_padic_tangent_p11.tsv",
    )
    ap.add_argument(
        "--p13",
        default="results/target_22224_direct_contact_deep13_padic_tangent_p13.tsv",
    )
    args = ap.parse_args()

    keys11 = deep_keys(Path(args.p11), 11, 11**5)
    keys13 = deep_keys(Path(args.p13), 13, 13**5)
    counts = {"rows": 0, "deep11": 0, "deep13": 0, "exact": 0, "full": 0}
    seen: set[tuple[int, ...]] = set()
    with Path(args.input).open() as src, Path(args.output).open("w") as dst:
        reader = csv.DictReader(src, delimiter="\t")
        dst.write(
            "left_index\tright_index\tsigns\tq1_num\tq1_den\tq2_num\tq2_den\t"
            "q3_num\tq3_den\ta\tb\tc\td\tsquare_mask\tfull_cover\n"
        )
        for row in reader:
            counts["rows"] += 1
            q = [(int(row[f"q{i}_num"]), int(row[f"q{i}_den"])) for i in range(1, 4)]
            signs = int(row["signs"])
            D = math.lcm(*(d for _, d in q))
            v = primitive(
                (D,)
                + tuple(((-n if signs >> i & 1 else n) * (D // d)) for i, (n, d) in enumerate(q))
            )
            if v in seen:
                continue
            seen.add(v)
            if canonical(v, 11, 11**5) not in keys11:
                raise AssertionError(("deep11", row, v))
            counts["deep11"] += 1
            if canonical(v, 13, 13**5) not in keys13:
                raise AssertionError(("deep13", row, v))
            counts["deep13"] += 1
            rr = radicands(v)
            mask = sum((1 << i) for i, z in enumerate(rr) if square(z))
            counts["exact"] += 1
            full = mask == 15
            counts["full"] += full
            dst.write(
                "\t".join(
                    str(row[k])
                    for k in (
                        "left_index",
                        "right_index",
                        "signs",
                        "q1_num",
                        "q1_den",
                        "q2_num",
                        "q2_den",
                        "q3_num",
                        "q3_den",
                    )
                )
                + "\t"
                + "\t".join(map(str, (*v, mask, int(full))))
                + "\n"
            )
    print(
        "A2228_DEEP_P5_CRT_VERIFY_DONE",
        counts,
        "unique",
        len(seen),
        "p11_keys",
        len(keys11),
        "p13_keys",
        len(keys13),
        "output",
        args.output,
    )


if __name__ == "__main__":
    main()
