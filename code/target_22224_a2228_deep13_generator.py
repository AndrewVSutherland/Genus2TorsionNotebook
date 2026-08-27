#!/usr/bin/env python3
"""Resumable rational-point generator for the corrected A2228+3 search.

The original version of this search used a faulty cubic-contact mask and
therefore imposed spurious ``deep p=13'' congruences.  The robust corrected
necessary condition is instead:

* bad branch reduction at both 11 and 13 (the corrected smooth intersections
  are empty); and
* at every other profiled prime, either bad branch reduction or
  3 | #J(F_p).

This script searches three genuine rational curves on the *full*
four-radicand A(2,2,2,8) cover: Filip 1, Filip 2, and Adam.  It consumes the
independently computed point-count profile TSV written by
target_22224_a2228_deep13_generator_profiles.m.  Work is chunked by the
denominator n, so disjoint jobs and exact resumption are straightforward.

Although the parameter height is modest, the homogeneous tuple coordinates
have degree four or six.  Thus even a height-100 parameter run already goes
far beyond the old integer bank bound B=16384 in tuple height.
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from math import gcd, isqrt
from pathlib import Path


FAMILY_NAMES = {1: "filip1", 2: "filip2", 3: "adam"}


def family_tuple(family: int, m: int, n: int) -> tuple[int, int, int, int]:
    """Homogeneous integral representatives for t=m/n."""
    if family == 1:
        return (
            -(m * m + m * n + n * n) ** 2,
            -4 * m * n * (m + n) ** 2,
            4 * m * n * n * (m + n),
            4 * m * m * n * (m + n),
        )
    if family == 2:
        h = m**4 + 2 * m**3 * n - m * m * n * n - 2 * m * n**3 + n**4
        z = m**4 - 2 * m**3 * n - m * m * n * n + 2 * m * n**3 + n**4
        return (-z * m * n, -h * n * n, h * m * n, h * m * m)
    if family == 3:
        return (
            -m * n * (m + n) ** 2,
            n * n * (m * m - n * n),
            -m * m * (m * m - n * n),
            m * n * (m - n) ** 2,
        )
    raise ValueError(f"unknown family {family}")


def primitive_tuple(vals: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for z in vals:
        g = gcd(g, abs(z))
    if g == 0:
        return vals
    out = tuple(z // g for z in vals)
    first = next((z for z in out if z), 1)
    return tuple(-z for z in out) if first < 0 else out


def cover_radicands(vals: tuple[int, int, int, int]) -> tuple[int, ...]:
    a, b, c, d = vals
    return (
        a * b * c * d,
        a * (a + b) * (a + c) * (a + d),
        b * (b + a) * (b + c) * (b + d),
        c * (c + a) * (c + b) * (c + d),
    )


def square_int(z: int) -> bool:
    if z < 0:
        return False
    r = isqrt(z)
    return r * r == z


def smooth_tuple(vals: tuple[int, int, int, int]) -> bool:
    return all(vals) and len({0, *(z * z for z in vals)}) == 5


def boundary_mod(vals: tuple[int, int, int, int], p: int) -> bool:
    return len({0, *((z * z) % p for z in vals)}) < 5


def curve_key(vals: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    sq = [z * z for z in vals]
    g = 0
    for z in sq:
        g = gcd(g, z)
    return tuple(sorted(z // g for z in sq))


@dataclass
class Profile:
    primes: tuple[int, ...]
    tables: dict[tuple[int, int], list[list[bool]]]
    allowed_counts: dict[tuple[int, int], int]


def load_profile(path: Path, requested_primes: set[int] | None) -> Profile:
    allowed_params: dict[tuple[int, int], set[int | str]] = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            family = int(row["family"])
            p = int(row["p"])
            if requested_primes is not None and p not in requested_primes:
                continue
            key = (family, p)
            allowed_params.setdefault(key, set())
            if int(row["allowed"]):
                label = row["parameter"]
                allowed_params[key].add("inf" if label == "inf" else int(label))

    tables: dict[tuple[int, int], list[list[bool]]] = {}
    counts: dict[tuple[int, int], int] = {}
    all_primes = sorted({p for _, p in allowed_params})
    for (family, p), params in allowed_params.items():
        table = [[False] * p for _ in range(p)]
        for mr in range(p):
            for nr in range(p):
                if mr == 0 and nr == 0:
                    continue
                if nr == 0:
                    ok = "inf" in params
                else:
                    r = mr * pow(nr, -1, p) % p
                    ok = r in params
                table[mr][nr] = ok
        tables[(family, p)] = table
        counts[(family, p)] = len(params)
    return Profile(tuple(all_primes), tables, counts)


def parse_int_set(text: str | None) -> set[int] | None:
    if text is None or text.strip().lower() in {"", "all"}:
        return None
    return {int(z) for z in text.split(",")}


def scan(args: argparse.Namespace) -> dict[str, int]:
    requested = parse_int_set(args.primes)
    profile = load_profile(Path(args.profiles), requested)
    families = sorted(parse_int_set(args.families) or set(FAMILY_NAMES))
    for family in families:
        missing = [p for p in profile.primes if (family, p) not in profile.tables]
        if missing:
            raise ValueError(f"family {family} missing profiles at {missing}")

    den_start = max(1, args.den_start)
    den_stop = args.den_stop if args.den_stop is not None else args.height
    checkpoint_path = Path(args.checkpoint) if args.checkpoint else None
    if args.resume and checkpoint_path and checkpoint_path.exists():
        saved = json.loads(checkpoint_path.read_text())
        den_start = max(den_start, int(saved["last_denominator"]) + 1)

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    append = args.resume and output.exists()
    handle = output.open("a" if append else "w", newline="")
    fields = [
        "family", "family_name", "m", "n", "a", "b", "c", "d",
        "tuple_height", "curve_key", "bad11", "bad13",
    ]
    writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields)
    if not append:
        writer.writeheader()

    tested = modular = smooth = cover = survivors = 0
    max_tuple_height = 0
    seen_curves: set[tuple[int, int, int, int]] = set()

    # Rarest parameter profiles first, which makes the inner loop cheap.
    prime_order: dict[int, list[int]] = {}
    for family in families:
        prime_order[family] = sorted(
            profile.primes,
            key=lambda p: profile.allowed_counts[(family, p)] / (p + 1),
        )
        print(
            "PROFILE_ORDER",
            FAMILY_NAMES[family],
            [(p, profile.allowed_counts[(family, p)]) for p in prime_order[family]],
            flush=True,
        )

    if den_start > den_stop:
        print("EMPTY_CHUNK", den_start, den_stop, flush=True)
    for n in range(den_start, den_stop + 1):
        nmods = {p: n % p for p in profile.primes}
        for family in families:
            ordered = prime_order[family]
            for m in range(-args.height, args.height + 1):
                tested += 1
                ok = True
                for p in ordered:
                    if not profile.tables[(family, p)][m % p][nmods[p]]:
                        ok = False
                        break
                if not ok or gcd(abs(m), n) != 1:
                    continue
                modular += 1
                vals = primitive_tuple(family_tuple(family, m, n))
                if not smooth_tuple(vals):
                    continue
                smooth += 1
                rads = cover_radicands(vals)
                if not all(square_int(z) for z in rads):
                    raise AssertionError((family, m, n, vals, rads))
                cover += 1
                if not boundary_mod(vals, 11) or not boundary_mod(vals, 13):
                    raise AssertionError(("missing forced boundary", family, m, n, vals))
                key = curve_key(vals)
                if key in seen_curves:
                    continue
                seen_curves.add(key)
                h = max(abs(z) for z in vals)
                max_tuple_height = max(max_tuple_height, h)
                writer.writerow(
                    {
                        "family": family,
                        "family_name": FAMILY_NAMES[family],
                        "m": m,
                        "n": n,
                        "a": vals[0],
                        "b": vals[1],
                        "c": vals[2],
                        "d": vals[3],
                        "tuple_height": h,
                        "curve_key": ",".join(map(str, key)),
                        "bad11": 1,
                        "bad13": 1,
                    }
                )
                handle.flush()
                survivors += 1
                print(
                    "MODULAR_SURVIVOR",
                    FAMILY_NAMES[family],
                    m,
                    n,
                    vals,
                    "height",
                    h,
                    flush=True,
                )
                if args.max_survivors and survivors >= args.max_survivors:
                    break
            if args.max_survivors and survivors >= args.max_survivors:
                break
        if checkpoint_path:
            checkpoint_path.write_text(
                json.dumps(
                    {
                        "last_denominator": n,
                        "height": args.height,
                        "profiles": str(args.profiles),
                        "primes": list(profile.primes),
                    },
                    indent=2,
                )
                + "\n"
            )
        if args.progress and n % args.progress == 0:
            print(
                "PROGRESS",
                n,
                "tested",
                tested,
                "modular",
                modular,
                "survivors",
                survivors,
                flush=True,
            )
        if args.max_survivors and survivors >= args.max_survivors:
            break

    handle.close()
    summary = {
        "tested": tested,
        "modular": modular,
        "smooth": smooth,
        "cover": cover,
        "survivors": survivors,
        "max_tuple_height": max_tuple_height,
        "den_start": den_start,
        "den_stop": den_stop,
    }
    print("SEARCH_DONE", json.dumps(summary, sort_keys=True), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--height", type=int, default=1000)
    parser.add_argument("--den-start", type=int, default=1)
    parser.add_argument("--den-stop", type=int)
    parser.add_argument("--families", default="1,2,3")
    parser.add_argument("--primes", default="all")
    parser.add_argument(
        "--profiles",
        default="results/target_22224_a2228_deep13_generator_profiles.tsv",
    )
    parser.add_argument(
        "--output",
        default="results/target_22224_a2228_deep13_generator_candidates.tsv",
    )
    parser.add_argument(
        "--checkpoint",
        default="results/target_22224_a2228_deep13_generator_checkpoint.json",
    )
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--progress", type=int, default=100)
    parser.add_argument("--max-survivors", type=int, default=0)
    args = parser.parse_args()
    scan(args)


if __name__ == "__main__":
    main()
