#!/usr/bin/env python3
"""Projective block-incidence masks for selected HLP pencil directions.

For a smooth genus-two fiber in odd characteristic, ``[6,12]`` implies the
three conditions below:

* rational 2-rank at least two;
* at least two distinct projective cubic-contact supports, equivalently
  rational 3-rank at least two;
* a rational half of a nonzero rational 2-class.

The contact and halving enumerators are the projective-chart enumerators from
``m612_hlp_slice_finite_masks.py``.  Smooth degree-five models are first
moved to degree six by x=a+1/z.  These masks are rigorous *necessary block
masks*, but they can have false positives on degenerate auxiliary supports;
the final answer must therefore use ``AbelianGroup`` and the actual invariant
factors.  The 72-divisibility and 2-rank checks are additional cheap
necessary prefilters.
"""

from __future__ import annotations

import argparse
import runpy
from dataclasses import dataclass, field
from pathlib import Path


HERE = Path(__file__).resolve().parent
ORDER = runpy.run_path(str(HERE / "m612_hlp_direction_order_prefilter.py"))
BLOCKS = runpy.run_path(str(HERE / "m612_hlp_slice_finite_masks.py"))
SEED = ORDER["SEED"]
PRIMES = (5, 7, 11, 13, 17, 19, 23)


def add(a: list[int], b: list[int], p: int) -> list[int]:
    out = [0] * max(len(a), len(b))
    for i in range(len(out)):
        out[i] = ((a[i] if i < len(a) else 0)
                  + (b[i] if i < len(b) else 0)) % p
    return ORDER["trim"](out, p)


def mul(a: list[int], b: list[int], p: int) -> list[int]:
    out = [0] * (len(a) + len(b) - 1)
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            out[i + j] = (out[i + j] + ai * bj) % p
    return ORDER["trim"](out, p)


def power(a: list[int], n: int, p: int) -> list[int]:
    result = [1]
    for _ in range(n):
        result = mul(result, a, p)
    return result


def degree_six_model(f: list[int], p: int) -> list[int]:
    """Return z^6 f(a+1/z), choosing f(a) nonzero, if deg(f)=5."""
    f = ORDER["trim"](f, p)
    if len(f) - 1 == 6:
        return f
    assert len(f) - 1 == 5
    a = next(x for x in range(p)
             if sum(c * pow(x, i, p) for i, c in enumerate(f)) % p)
    az1 = [1, a]  # 1+a*z
    result: list[int] = []
    for i, coefficient in enumerate(f):
        term = [0] * (6 - i) + [coefficient]
        term = mul(term, power(az1, i, p), p)
        result = add(result, term, p)
    assert len(result) - 1 == 6
    assert ORDER["is_smooth_genus2"](result, p)
    return result


def exact_status(f: list[int], p: int):
    """Return necessary block status: 1 passes, 0 killed, -1 bad, -2 error."""
    if not ORDER["is_smooth_genus2"](f, p):
        return -1, None
    order = ORDER["jacobian_order"](f, p)
    rank2 = ORDER["two_rank"](f, p)
    if order % 72 or rank2 < 2:
        return 0, None
    try:
        model = degree_six_model(f, p)
        contacts = BLOCKS["contact_supports"](model, p)
        if len(contacts) < 2:
            return 0, (order, rank2, len(contacts), 0)
        halves = BLOCKS["halving_supports"](model, p)
        record = (order, rank2, len(contacts), len(halves))
        return (1 if halves else 0), record
    except Exception as exc:  # retained separately; never silently killed
        return -2, repr(exc)


@dataclass
class Record:
    g: tuple[int, ...]
    allowed: list[tuple[int, ...]] = field(default_factory=list)
    bad: list[tuple[int, ...]] = field(default_factory=list)
    unresolved: list[tuple[int, ...]] = field(default_factory=list)
    details: list[tuple[tuple[int, object], ...]] = field(default_factory=list)

    def score(self):
        extras = [sum(t != 0 for t in mask) for mask in self.allowed]
        denominators = [
            (p - 1)
            - sum(t != 0 for t in bad)
            - sum(t != 0 for t in unresolved)
            for p, bad, unresolved in zip(PRIMES, self.bad, self.unresolved)
        ]
        breadth = sum(n > 0 for n in extras)
        density = sum(n / d for n, d in zip(extras, denominators) if d)
        return breadth, sum(extras), density, -sum(map(len, self.bad)), self.g


def read_directions(path: Path):
    seen = set()
    for raw in path.read_text().splitlines():
        raw = raw.strip()
        if not raw or raw.startswith("#"):
            continue
        g = tuple(map(int, raw.split(",")))
        if g not in seen:
            seen.add(g)
            yield g


def evaluate(record: Record) -> None:
    for p in PRIMES:
        allowed, bad, unresolved, details = [], [], [], []
        for t in range(p):
            f = [(SEED[i] + t * record.g[i]) % p for i in range(7)]
            status, detail = exact_status(f, p)
            if status == 1:
                allowed.append(t)
                details.append((t, detail))
            elif status == -1:
                bad.append(t)
            elif status == -2:
                unresolved.append(t)
                details.append((t, detail))
        record.allowed.append(tuple(allowed))
        record.bad.append(tuple(bad))
        record.unresolved.append(tuple(unresolved))
        record.details.append(tuple(details))


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--directions", type=Path, required=True)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--report-top", type=int, default=20)
    args = ap.parse_args()

    records = [Record(g) for g in read_directions(args.directions)]
    for i, record in enumerate(records, 1):
        evaluate(record)
        if i % 10 == 0:
            print(f"PROGRESS {i}/{len(records)}", flush=True)
    records.sort(key=Record.score, reverse=True)

    lines = [
        f"M612_HLP_DIRECTION_EXACT_BLOCK_SEARCH directions={len(records)}",
        "# infinity is retained separately and is not scored",
    ]
    for rank, record in enumerate(records[: args.report_top], 1):
        lines.append(f"RANK {rank} G={record.g} score={record.score()[:-1]}")
        for p, allowed, bad, unresolved, details in zip(
                PRIMES, record.allowed, record.bad, record.unresolved,
                record.details):
            lines.append(
                f" MASK p={p} allowed={allowed} bad={bad} "
                f"unresolved={unresolved} infinity=SEPARATE"
            )
            lines.append(f"  ALLOWED_RECORDS {details}")
    baseline = next((r for r in records if r.g == (1, 1, 0, 0, 0, 0, 0)), None)
    if baseline is not None:
        lines.append(f"BASELINE score={baseline.score()[:-1]}")
        for p, allowed, bad, unresolved in zip(
                PRIMES, baseline.allowed, baseline.bad, baseline.unresolved):
            lines.append(
                f" MASK p={p} allowed={allowed} bad={bad} "
                f"unresolved={unresolved} infinity=SEPARATE"
            )
    lines.append("M612_HLP_DIRECTION_EXACT_BLOCK_SEARCH_DONE")
    text = "\n".join(lines) + "\n"
    print(text, end="")
    if args.output:
        args.output.write_text(text)


if __name__ == "__main__":
    main()
