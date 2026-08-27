#!/usr/bin/env python3
"""p-adic residue analysis for the A(2,2,4,4) boundary square conditions.

This script works directly with residue classes on the K3 surface

    (ab + ac + ad + bc + bd + cd)^2 = 4abcd

and refines mod-p boundary ambiguities to mod-p^2.  For a fixed 2+2
partition, the A(2,2,4,4) local square conditions are tested by requiring
the four cross-differences a_i^2 - a_j^2 to have the same Q_p squareclass.
If a cross-difference is 0 modulo p^2, the class is recorded as a deep
collision rather than declared good.
"""

from __future__ import annotations

import argparse
import time
from collections import Counter
from itertools import product
from pathlib import Path


PARTITIONS = [
    ("12|34", (0, 1), (2, 3)),
    ("13|24", (0, 2), (1, 3)),
    ("14|23", (0, 3), (1, 2)),
]

PAIR_LABELS = {
    (0, 1): "E12",
    (0, 2): "E13",
    (0, 3): "E14",
    (1, 2): "E23",
    (1, 3): "E24",
    (2, 3): "E34",
}


def k3_F(x: tuple[int, int, int, int]) -> int:
    a, b, c, d = x
    s = a * b + a * c + a * d + b * c + b * d + c * d
    return s * s - 4 * a * b * c * d


def k3_grad(x: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    a, b, c, d = x
    s = a * b + a * c + a * d + b * c + b * d + c * d
    return (
        2 * s * (b + c + d) - 4 * b * c * d,
        2 * s * (a + c + d) - 4 * a * c * d,
        2 * s * (a + b + d) - 4 * a * b * d,
        2 * s * (a + b + c) - 4 * a * b * c,
    )


def legendre(u: int, p: int) -> int:
    u %= p
    if u == 0:
        raise ValueError("Legendre symbol requested for zero unit")
    return 1 if pow(u, (p - 1) // 2, p) == 1 else -1


def squareclass_table(p: int, depth: int) -> list[tuple[int, int] | None]:
    q = p**depth
    table: list[tuple[int, int] | None] = [None] * q
    for r in range(1, q):
        n = r
        v = 0
        while v < depth and n % p == 0:
            v += 1
            n //= p
        if v == depth:
            table[r] = None
        else:
            table[r] = (v % 2, legendre(n, p))
    return table


def cross_pairs(left: tuple[int, int], right: tuple[int, int]) -> list[tuple[int, int]]:
    return [(i, j) for i in left for j in right]


def classify_partition(
    squares: list[int],
    q: int,
    sc_table: list[tuple[int, int] | None],
    pairs: list[tuple[int, int]],
) -> tuple[str, tuple[str, ...]]:
    known: set[tuple[int, int]] = set()
    deep: list[str] = []
    for i, j in pairs:
        sc = sc_table[(squares[i] - squares[j]) % q]
        if sc is None:
            deep.append(PAIR_LABELS[(i, j) if i < j else (j, i)])
        else:
            known.add(sc)
    if len(known) > 1:
        return "killed", ()
    if deep:
        return "deep", tuple(deep)
    return "resolved_ok", ()


def boundary_signature(x: tuple[int, int, int, int], p: int) -> tuple[str, ...]:
    labels: list[str] = []
    squares = [(v * v) % p for v in x]
    for i, v in enumerate(x):
        if v % p == 0:
            labels.append(f"Z{i + 1}")
    for i in range(4):
        for j in range(i + 1, 4):
            if squares[i] == squares[j]:
                labels.append(PAIR_LABELS[(i, j)])
    return tuple(labels)


def lift_vectors_mod_p2(x0: tuple[int, int, int, int], p: int):
    f0 = k3_F(x0)
    if f0 % p != 0:
        raise ValueError("base point is not on the K3 surface mod p")

    rhs = (-(f0 // p)) % p
    grad = [g % p for g in k3_grad(x0)]
    pivot = next((i for i, g in enumerate(grad) if g), None)

    if pivot is None:
        if f0 % (p * p) != 0:
            return
        for xs in product(range(p), repeat=4):
            yield xs
        return

    inv = pow(grad[pivot], -1, p)
    free = [i for i in range(4) if i != pivot]
    for vals in product(range(p), repeat=3):
        lift = [0, 0, 0, 0]
        total = 0
        for i, val in zip(free, vals):
            lift[i] = val
            total += grad[i] * val
        lift[pivot] = ((rhs - total) * inv) % p
        yield tuple(lift)


def analyze_prime(p: int, max_examples: int, include_origin: bool) -> dict:
    t0 = time.time()
    q = p * p
    sc_p = squareclass_table(p, 1)
    sc_p2 = squareclass_table(p, 2)
    partition_pairs = {
        name: cross_pairs(left, right) for name, left, right in PARTITIONS
    }

    summary: dict = {
        "p": p,
        "roots": 0,
        "boundary_roots": 0,
        "good_roots": 0,
        "singular_roots": 0,
        "good_roots_with_any_modp_survivor": 0,
        "ambiguous_roots": 0,
        "modp_partition_status": {name: Counter() for name, _, _ in PARTITIONS},
        "modp_root_status_distribution": Counter(),
        "ambiguous_boundary_signatures": Counter(),
        "p2_lift_status_by_partition": {name: Counter() for name, _, _ in PARTITIONS},
        "p2_root_partition_outcomes": {name: Counter() for name, _, _ in PARTITIONS},
        "p2_deep_masks": {name: Counter() for name, _, _ in PARTITIONS},
        "p2_resolved_boundary_signatures": {name: Counter() for name, _, _ in PARTITIONS},
        "p2_deep_boundary_signatures": {name: Counter() for name, _, _ in PARTITIONS},
        "examples": {"resolved_ok": [], "deep": [], "all_lifts_killed": []},
    }

    ambiguous: list[tuple[tuple[int, int, int, int], list[str], tuple[str, ...], bool]] = []

    for x0 in product(range(p), repeat=4):
        if not include_origin and all(v == 0 for v in x0):
            continue
        if k3_F(x0) % p != 0:
            continue

        summary["roots"] += 1
        grad = tuple(g % p for g in k3_grad(x0))
        singular = all(g == 0 for g in grad)
        if singular:
            summary["singular_roots"] += 1

        sig = boundary_signature(x0, p)
        if sig:
            summary["boundary_roots"] += 1
        else:
            summary["good_roots"] += 1

        squares = [(v * v) % p for v in x0]
        statuses: list[tuple[str, str]] = []
        ambiguous_parts: list[str] = []
        has_survivor = False
        for name, _, _ in PARTITIONS:
            status, _ = classify_partition(squares, p, sc_p, partition_pairs[name])
            if status == "deep":
                status = "ambiguous"
                ambiguous_parts.append(name)
            if status == "resolved_ok":
                has_survivor = True
            summary["modp_partition_status"][name][status] += 1
            statuses.append((name, status))

        if has_survivor and not sig:
            summary["good_roots_with_any_modp_survivor"] += 1
        summary["modp_root_status_distribution"][
            tuple(status for _, status in statuses)
        ] += 1
        if ambiguous_parts:
            summary["ambiguous_roots"] += 1
            for name in ambiguous_parts:
                summary["ambiguous_boundary_signatures"][(name, sig)] += 1
            ambiguous.append((x0, ambiguous_parts, sig, singular))

    summary["ambiguous_root_partitions"] = sum(
        c["ambiguous"] for c in summary["modp_partition_status"].values()
    )

    for root_index, (x0, ambiguous_parts, sig, singular) in enumerate(ambiguous, 1):
        part_counts = {name: Counter() for name in ambiguous_parts}
        part_deep_masks = {name: Counter() for name in ambiguous_parts}
        first_example: dict[str, tuple[int, int, int, int] | None] = {
            name: None for name in ambiguous_parts
        }
        first_deep: dict[str, tuple[tuple[int, int, int, int], tuple[str, ...]] | None] = {
            name: None for name in ambiguous_parts
        }

        for lift in lift_vectors_mod_p2(x0, p):
            x = tuple(x0[i] + p * lift[i] for i in range(4))
            squares = [(v * v) % q for v in x]
            for name in ambiguous_parts:
                status, deep_mask = classify_partition(
                    squares, q, sc_p2, partition_pairs[name]
                )
                summary["p2_lift_status_by_partition"][name][status] += 1
                part_counts[name][status] += 1
                if status == "deep":
                    part_deep_masks[name][deep_mask] += 1
                    summary["p2_deep_masks"][name][deep_mask] += 1
                    summary["p2_deep_boundary_signatures"][name][sig] += 1
                    if first_deep[name] is None:
                        first_deep[name] = (x, deep_mask)
                elif status == "resolved_ok":
                    summary["p2_resolved_boundary_signatures"][name][sig] += 1
                    if first_example[name] is None:
                        first_example[name] = x

        for name in ambiguous_parts:
            counts = part_counts[name]
            if counts["resolved_ok"] and counts["deep"]:
                outcome = "resolved_and_deep"
            elif counts["resolved_ok"]:
                outcome = "resolved_only"
            elif counts["deep"]:
                outcome = "deep_only"
            else:
                outcome = "all_lifts_killed"
            summary["p2_root_partition_outcomes"][name][outcome] += 1
            if counts["resolved_ok"]:
                summary["p2_root_partition_outcomes"][name]["has_resolved_ok_lift"] += 1
            if counts["deep"]:
                summary["p2_root_partition_outcomes"][name]["has_deep_lift"] += 1

            if (
                outcome == "all_lifts_killed"
                and len(summary["examples"]["all_lifts_killed"]) < max_examples
            ):
                summary["examples"]["all_lifts_killed"].append((p, x0, name, sig, singular))
            if (
                first_example[name] is not None
                and len(summary["examples"]["resolved_ok"]) < max_examples
            ):
                summary["examples"]["resolved_ok"].append(
                    (p, x0, first_example[name], name, sig, singular)
                )
            if (
                first_deep[name] is not None
                and len(summary["examples"]["deep"]) < max_examples
            ):
                x, deep_mask = first_deep[name]
                summary["examples"]["deep"].append(
                    (p, x0, x, name, sig, deep_mask, singular)
                )

        if root_index % 500 == 0:
            elapsed = time.time() - t0
            print(
                f"p={p}: refined {root_index}/{len(ambiguous)} ambiguous roots "
                f"({elapsed:.1f}s)",
                flush=True,
            )

    summary["elapsed_seconds"] = time.time() - t0
    return summary


def format_counter(counter: Counter, keys: list[str] | None = None) -> str:
    if keys is None:
        keys = sorted(counter)
    return " ".join(f"{key}={counter[key]}" for key in keys)


def write_report(path: Path, summaries: list[dict], top: int) -> None:
    with path.open("w") as out:
        out.write("A2244 p-adic residue-class boundary analysis\n")
        out.write("surface=(ab+ac+ad+bc+bd+cd)^2-4abcd=0\n")
        out.write("depth=p^2 refinement of mod-p ambiguous partition classes\n")
        out.write("all_zero_mod_p_root=excluded_by_default\n")
        out.write("\nStatus meanings\n")
        out.write("resolved_ok: all four cross-differences have a determined common squareclass\n")
        out.write("deep: compatible so far, but at least one cross-difference is 0 mod p^2\n")
        out.write("killed: determined cross-difference squareclasses are incompatible\n")

        for summary in summaries:
            p = summary["p"]
            out.write(f"\np={p}\n")
            out.write(f"elapsed_seconds={summary['elapsed_seconds']:.3f}\n")
            out.write(f"Fp_roots={summary['roots']}\n")
            out.write(f"boundary_roots={summary['boundary_roots']}\n")
            out.write(f"good_reduction_roots={summary['good_roots']}\n")
            out.write(f"singular_roots={summary['singular_roots']}\n")
            out.write(f"ambiguous_roots={summary['ambiguous_roots']}\n")
            out.write(f"ambiguous_root_partitions={summary['ambiguous_root_partitions']}\n")
            out.write(
                "good_reduction_roots_with_any_modp_partition_survivor="
                f"{summary['good_roots_with_any_modp_survivor']}\n"
            )

            out.write("\nmod-p partition status\n")
            for name, _, _ in PARTITIONS:
                out.write(
                    f"  {name} "
                    + format_counter(
                        summary["modp_partition_status"][name],
                        ["killed", "ambiguous", "resolved_ok"],
                    )
                    + "\n"
                )

            out.write("\nmod-p root status distributions\n")
            for statuses, count in summary["modp_root_status_distribution"].most_common():
                out.write(f"  {count} {statuses}\n")

            out.write("\np^2 lift residue status, restricted to mod-p ambiguous root-partitions\n")
            for name, _, _ in PARTITIONS:
                out.write(
                    f"  {name} "
                    + format_counter(
                        summary["p2_lift_status_by_partition"][name],
                        ["killed", "deep", "resolved_ok"],
                    )
                    + "\n"
                )

            out.write("\np^2 root-partition outcomes\n")
            for name, _, _ in PARTITIONS:
                out.write(
                    f"  {name} "
                    + format_counter(
                        summary["p2_root_partition_outcomes"][name],
                        [
                            "all_lifts_killed",
                            "deep_only",
                            "resolved_only",
                            "resolved_and_deep",
                            "has_resolved_ok_lift",
                            "has_deep_lift",
                        ],
                    )
                    + "\n"
                )

            out.write("\ntop mod-p boundary signatures among ambiguous root-partitions\n")
            for (name, sig), count in summary["ambiguous_boundary_signatures"].most_common(top):
                out.write(f"  {count} {name} {'+'.join(sig) if sig else 'good'}\n")

            out.write("\ntop p^2 deep masks by partition\n")
            for name, _, _ in PARTITIONS:
                out.write(f"  {name}\n")
                for mask, count in summary["p2_deep_masks"][name].most_common(top):
                    out.write(f"    {count} {'+'.join(mask)}\n")

            out.write("\ntop mod-p boundary signatures with p^2 resolved lifts\n")
            for name, _, _ in PARTITIONS:
                out.write(f"  {name}\n")
                for sig, count in summary["p2_resolved_boundary_signatures"][name].most_common(top):
                    out.write(f"    {count} {'+'.join(sig) if sig else 'good'}\n")

            out.write("\ntop mod-p boundary signatures with p^2 deep lifts\n")
            for name, _, _ in PARTITIONS:
                out.write(f"  {name}\n")
                for sig, count in summary["p2_deep_boundary_signatures"][name].most_common(top):
                    out.write(f"    {count} {'+'.join(sig) if sig else 'good'}\n")

            out.write("\nexamples: resolved_ok p^2 lifts\n")
            for item in summary["examples"]["resolved_ok"]:
                _, x0, x, name, sig, singular = item
                out.write(
                    f"  modp={x0} lift={x} partition={name} "
                    f"boundary={' + '.join(sig) if sig else 'good'} singular={singular}\n"
                )

            out.write("\nexamples: deep p^2 lifts\n")
            for item in summary["examples"]["deep"]:
                _, x0, x, name, sig, deep_mask, singular = item
                out.write(
                    f"  modp={x0} lift={x} partition={name} "
                    f"deep={'+'.join(deep_mask)} "
                    f"boundary={' + '.join(sig) if sig else 'good'} singular={singular}\n"
                )

            out.write("\nexamples: mod-p ambiguous root-partitions killed by every p^2 lift\n")
            for item in summary["examples"]["all_lifts_killed"]:
                _, x0, name, sig, singular = item
                out.write(
                    f"  modp={x0} partition={name} "
                    f"boundary={' + '.join(sig) if sig else 'good'} singular={singular}\n"
                )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", nargs="+", type=int, default=[11, 23])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/a2244_padic_residue_report.txt"),
    )
    parser.add_argument("--top", type=int, default=12)
    parser.add_argument("--max-examples", type=int, default=12)
    parser.add_argument(
        "--include-origin",
        action="store_true",
        help="include the all-zero mod-p root; off by default for primitive projective tuples",
    )
    args = parser.parse_args()

    summaries = []
    for p in args.primes:
        print(f"analyzing p={p}", flush=True)
        summaries.append(analyze_prime(p, args.max_examples, args.include_origin))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_report(args.output, summaries, args.top)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
