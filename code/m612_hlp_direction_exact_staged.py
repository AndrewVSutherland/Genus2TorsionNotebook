#!/usr/bin/env python3
"""Deterministic staged projective-block search over the small direction box.

The necessary block-incidence fiber predicate is imported from
``m612_hlp_direction_exact_block_search.py``.  Projectively congruent
directions are cached at each prime.  Staging keeps the search bounded:
all directions are tested at 7, the best ``stage_top`` continue to 11 and
13, and the best ``final_top`` continue to 17, 19, and 23.  Prime 5 is
then evaluated for the finalists solely to retain its bad-fiber mask; no
box-two direction has a smooth 72-divisible fiber there.

The output is a finalist generator, not the final finite-group answer:
leading directions must be checked by ``m612_hlp_direction_search.m``.
"""

from __future__ import annotations

import argparse
import runpy
from dataclasses import dataclass, field
from pathlib import Path


HERE = Path(__file__).resolve().parent
EXACT = runpy.run_path(str(HERE / "m612_hlp_direction_exact_block_search.py"))
ORDER = EXACT["ORDER"]
SEED = EXACT["SEED"]
PRIMES = (5, 7, 11, 13, 17, 19, 23)


@dataclass
class Record:
    g: tuple[int, ...]
    allowed: dict[int, tuple[int, ...]] = field(default_factory=dict)
    bad: dict[int, tuple[int, ...]] = field(default_factory=dict)
    unresolved: dict[int, tuple[int, ...]] = field(default_factory=dict)
    details: dict[int, tuple[tuple[int, object], ...]] = field(default_factory=dict)

    def score(self, primes: tuple[int, ...]):
        extras = [sum(t != 0 for t in self.allowed[p]) for p in primes]
        denominators = [
            (p - 1)
            - sum(t != 0 for t in self.bad[p])
            - sum(t != 0 for t in self.unresolved[p])
            for p in primes
        ]
        breadth = sum(n > 0 for n in extras)
        density = sum(n / d for n, d in zip(extras, denominators) if d)
        return breadth, sum(extras), density, -sum(len(self.bad[p]) for p in primes), self.g


def normalized(g: tuple[int, ...], p: int):
    w = tuple(a % p for a in g)
    scale = next(a for a in w if a)
    inverse = pow(scale, -1, p)
    return tuple(a * inverse % p for a in w), scale


def base_line_mask(n: tuple[int, ...], p: int, seed_info):
    allowed, bad, unresolved, details = [], [], [], []
    status, detail = seed_info
    if status == 1:
        allowed.append(0)
        details.append((0, detail))
    elif status == -1:
        bad.append(0)
    elif status == -2:
        unresolved.append(0)
        details.append((0, detail))
    for u in range(1, p):
        f = [(SEED[i] + u * n[i]) % p for i in range(7)]
        status, detail = EXACT["exact_status"](f, p)
        if status == 1:
            allowed.append(u)
            details.append((u, detail))
        elif status == -1:
            bad.append(u)
        elif status == -2:
            unresolved.append(u)
            details.append((u, detail))
    return tuple(allowed), tuple(bad), tuple(unresolved), tuple(details)


def add_prime(records: list[Record], p: int) -> None:
    cache = {}
    seed_f = [a % p for a in SEED]
    seed_info = EXACT["exact_status"](seed_f, p)
    for j, record in enumerate(records, 1):
        n, scale = normalized(record.g, p)
        if n not in cache:
            cache[n] = base_line_mask(n, p, seed_info)
        allowed, bad, unresolved, details = cache[n]
        inverse = pow(scale, -1, p)
        record.allowed[p] = tuple(sorted(u * inverse % p for u in allowed))
        record.bad[p] = tuple(sorted(u * inverse % p for u in bad))
        record.unresolved[p] = tuple(sorted(u * inverse % p for u in unresolved))
        record.details[p] = tuple(sorted(
            ((u * inverse % p, detail) for u, detail in details),
            key=lambda item: item[0],
        ))
        if j % 1000 == 0:
            print(f"PRIME_PROGRESS p={p} records={j}/{len(records)} "
                  f"cache={len(cache)}", flush=True)
    print(f"PRIME_DONE p={p} records={len(records)} "
          f"unique_projective={len(cache)}", flush=True)


def report(label: str, records: list[Record], primes: tuple[int, ...], n: int):
    print(f"REPORT {label} records={len(records)} showing={min(n, len(records))}")
    for rank, record in enumerate(records[:n], 1):
        print(f"RANK {rank} G={record.g} score={record.score(primes)[:-1]}")
        for p in primes:
            print(f" MASK p={p} allowed={record.allowed[p]} bad={record.bad[p]} "
                  f"unresolved={record.unresolved[p]} infinity=SEPARATE")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--box", type=int, default=2)
    ap.add_argument("--stage-top", type=int, default=5000)
    ap.add_argument("--final-top", type=int, default=1000)
    ap.add_argument("--report-top", type=int, default=30)
    ap.add_argument("--output", type=Path)
    args = ap.parse_args()

    records = [Record(tuple(g)) for g in ORDER["primitive_transverse_directions"](args.box)]
    print(f"M612_HLP_DIRECTION_EXACT_STAGED box={args.box} directions={len(records)}")

    add_prime(records, 7)
    records.sort(key=lambda r: r.score((7,)), reverse=True)
    report("STAGE1", records, (7,), min(20, args.report_top))
    records = records[: args.stage_top]

    for p in (11, 13):
        add_prime(records, p)
    records.sort(key=lambda r: r.score((7, 11, 13)), reverse=True)
    report("STAGE2", records, (7, 11, 13), min(20, args.report_top))
    records = records[: args.final_top]

    for p in (17, 19, 23):
        add_prime(records, p)
    records.sort(key=lambda r: r.score((7, 11, 13, 17, 19, 23)), reverse=True)
    add_prime(records, 5)

    # Capture the final report so it is both printed and reproducibly saved.
    lines = [
        f"M612_HLP_DIRECTION_EXACT_STAGED_FINAL records={len(records)}",
        "# infinity is separate; score uses p=7,11,13,17,19,23",
    ]
    for rank, record in enumerate(records[: args.report_top], 1):
        lines.append(
            f"RANK {rank} G={record.g} "
            f"score={record.score((7,11,13,17,19,23))[:-1]}"
        )
        for p in PRIMES:
            lines.append(
                f" MASK p={p} allowed={record.allowed[p]} bad={record.bad[p]} "
                f"unresolved={record.unresolved[p]} infinity=SEPARATE"
            )
            lines.append(f"  ALLOWED_RECORDS {record.details[p]}")
    lines.append("M612_HLP_DIRECTION_EXACT_STAGED_DONE")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    if args.output:
        args.output.write_text(text)


if __name__ == "__main__":
    main()
