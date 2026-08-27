#!/usr/bin/env python3
"""Small-prime congruence diagnostics for A(2,2,4,4).

For p <= 7 every primitive K3 residue class has bad reduction in the
four squared branch points: there are fewer than four nonzero square
residues in F_p.  This script checks whether the finer p-adic
A(2,2,4,4) squareclass condition still gives useful congruence filters.

The surface is

    (ab + ac + ad + bc + bd + cd)^2 = 4abcd.

For a fixed 2+2 partition, the A(2,2,4,4) condition is that the four
oriented cross-differences a_i^2 - a_j^2 have the same Q_p squareclass.
For odd p this is valuation parity plus Legendre symbol of the unit.  For
p = 2 it is valuation parity plus odd unit modulo 8.
"""

from __future__ import annotations

import argparse
from collections import Counter
from itertools import product
from pathlib import Path


PARTITIONS = [
    ("12|34", (0, 1), (2, 3)),
    ("13|24", (0, 2), (1, 3)),
    ("14|23", (0, 3), (1, 2)),
]

PAIR_LABELS = ["E12", "E13", "E14", "E23", "E24", "E34"]


def k3_F(x: tuple[int, int, int, int]) -> int:
    a, b, c, d = x
    s = a * b + a * c + a * d + b * c + b * d + c * d
    return s * s - 4 * a * b * c * d


def legendre(u: int, p: int) -> int:
    u %= p
    if u == 0:
        raise ValueError("Legendre symbol requested for zero")
    return 1 if pow(u, (p - 1) // 2, p) == 1 else -1


def all_q2_unit_classes(parity: int) -> set[tuple[int, int]]:
    return {(parity, u) for u in (1, 3, 5, 7)}


def squareclass_possibilities(n: int, p: int, depth: int) -> set[tuple[int, int]] | None:
    """Return possible Q_p squareclasses from n modulo p^depth.

    None means n is 0 modulo p^depth, so the finite residue data are still
    deep and impose no determined squareclass at this precision.
    """

    q = p**depth
    r = n % q
    if r == 0:
        return None

    v = 0
    while r % p == 0:
        v += 1
        r //= p

    parity = v & 1
    if p != 2:
        return {(parity, legendre(r, p))}

    known_unit_bits = depth - v
    if known_unit_bits >= 3:
        return {(parity, r % 8)}
    if known_unit_bits == 2:
        return {(parity, u) for u in (1, 3, 5, 7) if u % 4 == r % 4}
    return all_q2_unit_classes(parity)


def partition_status(
    t: tuple[int, int, int, int],
    p: int,
    depth: int,
    left: tuple[int, int],
    right: tuple[int, int],
) -> str:
    q = p**depth
    squares = [(x * x) % q for x in t]
    possible_common: set[tuple[int, int]] | None = None
    saw_deep = False
    saw_ambiguous = False

    for i in left:
        for j in right:
            classes = squareclass_possibilities(squares[i] - squares[j], p, depth)
            if classes is None:
                saw_deep = True
                continue
            if len(classes) > 1:
                saw_ambiguous = True
            possible_common = classes if possible_common is None else possible_common & classes
            if not possible_common:
                return "killed"

    if saw_deep:
        return "deep"
    if saw_ambiguous:
        return "unit_ambiguous"
    return "resolved_ok"


def boundary_signature_mod_p(t: tuple[int, int, int, int], p: int) -> str:
    squares = [(x * x) % p for x in t]
    labels: list[str] = []
    for i, x in enumerate(t):
        if x % p == 0:
            labels.append(f"Z{i + 1}")

    idx = 0
    for i in range(4):
        for j in range(i + 1, 4):
            if squares[i] == squares[j]:
                labels.append(PAIR_LABELS[idx])
            idx += 1
    return "+".join(labels) if labels else "good"


def analyze_case(p: int, depth: int, top: int) -> dict:
    q = p**depth
    summary = {
        "p": p,
        "depth": depth,
        "modulus": q,
        "roots": 0,
        "primitive_roots": 0,
        "boundary_roots_mod_p": 0,
        "any_partition_not_killed": 0,
        "real_partition_12_34_not_killed": 0,
        "partition_status": {name: Counter() for name, _, _ in PARTITIONS},
        "signature_counts": Counter(),
        "possible_signature_counts": Counter(),
        "fully_killed_signature_counts": Counter(),
        "top": top,
    }

    for t in product(range(q), repeat=4):
        if k3_F(t) % q != 0:
            continue
        summary["roots"] += 1
        if all(x % p == 0 for x in t):
            continue
        summary["primitive_roots"] += 1

        sig = boundary_signature_mod_p(t, p)
        summary["signature_counts"][sig] += 1
        if sig != "good":
            summary["boundary_roots_mod_p"] += 1

        statuses: list[tuple[str, str]] = []
        for name, left, right in PARTITIONS:
            status = partition_status(t, p, depth, left, right)
            summary["partition_status"][name][status] += 1
            statuses.append((name, status))

        if any(status != "killed" for _, status in statuses):
            summary["any_partition_not_killed"] += 1
            summary["possible_signature_counts"][sig] += 1
        else:
            summary["fully_killed_signature_counts"][sig] += 1

        if statuses[0][1] != "killed":
            summary["real_partition_12_34_not_killed"] += 1

    return summary


def parse_cases(raw: str) -> list[tuple[int, int]]:
    out = []
    for part in raw.split(","):
        if not part:
            continue
        p_raw, depth_raw = part.split(":", 1)
        out.append((int(p_raw), int(depth_raw)))
    return out


def write_report(path: Path, summaries: list[dict]) -> None:
    with path.open("w") as out:
        out.write("A2244 small-prime bad-reduction congruence diagnostics\n")
        out.write("surface=(ab+ac+ad+bc+bd+cd)^2-4abcd=0\n\n")
        out.write("Bad reduction condition for squared branch coordinates:\n")
        out.write("  p | a*b*c*d*prod_{i<j}(a_i^2-a_j^2).\n")
        out.write("For p <= 7 this is automatic for every primitive tuple, since\n")
        out.write("|(F_p^*)^2|=(p-1)/2 < 4.  Equivalently, modulo p at least\n")
        out.write("one coordinate is zero or two coordinates are congruent up to sign.\n\n")
        out.write("A2244 squareclass condition for a partition I|J:\n")
        out.write("  all a_i^2-a_j^2, i in I, j in J, have the same Q_p squareclass.\n")
        out.write("For odd p: same valuation parity and same Legendre symbol of the unit.\n")
        out.write("For p=2: same valuation parity and same odd unit modulo 8.\n")
        out.write("Deep means a cross-difference is 0 modulo the tested prime power.\n")
        out.write("Unit_ambiguous only occurs at p=2 when the tested modulus has not\n")
        out.write("determined the odd unit modulo 8.\n")

        for summary in summaries:
            p = summary["p"]
            depth = summary["depth"]
            primitive = summary["primitive_roots"]
            out.write(f"\np={p}, depth={depth}, modulus={summary['modulus']}\n")
            out.write(f"roots_mod_prime_power={summary['roots']}\n")
            out.write(f"primitive_roots={primitive}\n")
            out.write(f"boundary_roots_mod_p={summary['boundary_roots_mod_p']}\n")
            out.write(
                "any_partition_not_killed="
                f"{summary['any_partition_not_killed']}"
                f" ({summary['any_partition_not_killed'] / primitive:.6f})\n"
            )
            out.write(
                "real_partition_12|34_not_killed="
                f"{summary['real_partition_12_34_not_killed']}"
                f" ({summary['real_partition_12_34_not_killed'] / primitive:.6f})\n"
            )
            out.write("partition_status_counts\n")
            for name, _, _ in PARTITIONS:
                counts = summary["partition_status"][name]
                out.write(
                    f"  {name} "
                    + " ".join(f"{key}={counts[key]}" for key in sorted(counts))
                    + "\n"
                )
            out.write("top_mod_p_boundary_signatures\n")
            for sig, count in summary["signature_counts"].most_common(summary["top"]):
                out.write(f"  {count} {sig}\n")
            out.write("top_not_killed_signatures\n")
            for sig, count in summary["possible_signature_counts"].most_common(summary["top"]):
                out.write(f"  {count} {sig}\n")
            out.write("top_fully_killed_signatures\n")
            for sig, count in summary["fully_killed_signature_counts"].most_common(summary["top"]):
                out.write(f"  {count} {sig}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--cases",
        default="2:5,3:3,5:2,7:1",
        help="comma-separated p:depth cases",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/a2244_small_prime_congruences.txt"),
    )
    parser.add_argument("--top", type=int, default=10)
    args = parser.parse_args()

    summaries = []
    for p, depth in parse_cases(args.cases):
        print(f"analyzing p={p}, depth={depth}", flush=True)
        summaries.append(analyze_case(p, depth, args.top))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_report(args.output, summaries)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
