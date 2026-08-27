#!/usr/bin/env python3
"""Sampled p^3 refinement of normalized A2244 p^2-deep boundary charts.

The exact p^2-deep residue population is too large to lift exhaustively to
p^3 at p=23.  This script compresses p^2-deep classes by a broad normalized
chart key and stores several control representatives in each chart.  The
control representative key records the p^3 K3 lift obstruction and the p^2
coefficients of the deep A(2,2,4,4) cross-differences.  It then lifts those
representatives to p^3.

This is a chart diagnostic: it identifies persistent and mixed local chart
behavior, but it is not a full enumeration of every p^2 residue class.
"""

from __future__ import annotations

import argparse
import time
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

import analyze_a2244_padic_residues as base

PAIR_ORDER = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
STATUS_ORDER = ["killed", "deep", "resolved_ok"]


def vp_cap(x: int, p: int, depth: int) -> int:
    q = p**depth
    x %= q
    if x == 0:
        return depth
    v = 0
    while v < depth and x % p == 0:
        v += 1
        x //= p
    return v


def pair_relations_mod_p2(x: tuple[int, int, int, int], p: int) -> tuple[str, ...]:
    q = p * p
    rels = []
    for i, j in PAIR_ORDER:
        if (x[i] * x[i] - x[j] * x[j]) % q:
            rel = "."
        elif x[i] % p == 0 and x[j] % p == 0:
            rel = "Z"
        elif (x[i] - x[j]) % q == 0:
            rel = "+"
        elif (x[i] + x[j]) % q == 0:
            rel = "-"
        else:
            rel = "?"
        rels.append(rel)
    return tuple(rels)


def known_cross_classes(squares, q, sc_table, pairs):
    known = set()
    for i, j in pairs:
        sc = sc_table[(squares[i] - squares[j]) % q]
        if sc is not None:
            known.add(sc)
    return tuple(sorted(known))




def classify_partition_with_known(squares, q, sc_table, pairs):
    known = set()
    deep = []
    for i, j in pairs:
        sc = sc_table[(squares[i] - squares[j]) % q]
        if sc is None:
            deep.append(base.PAIR_LABELS[(i, j) if i < j else (j, i)])
        else:
            known.add(sc)
    if len(known) > 1:
        return "killed", (), tuple(sorted(known))
    if deep:
        return "deep", tuple(deep), tuple(sorted(known))
    return "resolved_ok", (), tuple(sorted(known))


def deep_coeffs_mod_p2(x, p, pairs, deep_mask):
    q = p * p
    coeffs = []
    deep = set(deep_mask)
    for i, j in pairs:
        label = base.PAIR_LABELS[(i, j) if i < j else (j, i)]
        if label in deep:
            coeffs.append((label, ((x[i] * x[i] - x[j] * x[j]) // q) % p))
    return tuple(coeffs)


def coeff_signature(coeffs, p):
    sig = []
    for label, coeff in coeffs:
        if coeff == 0:
            sig.append((label, 0))
        else:
            sig.append((label, base.legendre(coeff, p)))
    return tuple(sig)


def lift_vectors_mod_p3_from_p2(x2, p):
    q2 = p * p
    q3 = q2 * p
    f2 = base.k3_F(x2)
    if f2 % q2 != 0:
        raise ValueError("p^2 representative is not on the K3 surface mod p^2")

    rhs = (-(f2 // q2)) % p
    grad = [g % p for g in base.k3_grad(x2)]
    pivot = next((i for i, g in enumerate(grad) if g), None)

    if pivot is None:
        if f2 % q3 != 0:
            return
        for vals in product(range(p), repeat=4):
            yield vals
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


def modp_ambiguous_roots(p):
    sc_p = base.squareclass_table(p, 1)
    pairs = {name: base.cross_pairs(left, right) for name, left, right in base.PARTITIONS}
    ambiguous = []
    for x0 in product(range(p), repeat=4):
        if all(v == 0 for v in x0):
            continue
        if base.k3_F(x0) % p:
            continue
        squares = [(v * v) % p for v in x0]
        parts = []
        for name, _, _ in base.PARTITIONS:
            status, _ = base.classify_partition(squares, p, sc_p, pairs[name])
            if status == "deep":
                parts.append(name)
        if parts:
            singular = all(g % p == 0 for g in base.k3_grad(x0))
            ambiguous.append((x0, parts, base.boundary_signature(x0, p), singular))
    return ambiguous


def broad_chart_key(p, partition, sig, deep_mask, x2, known, singular):
    return (
        partition,
        sig,
        deep_mask,
        tuple(vp_cap(v, p, 2) for v in x2),
        pair_relations_mod_p2(x2, p),
        known,
        singular,
    )


def control_key(p, x2, pairs, deep_mask):
    q2 = p * p
    grad = tuple(g % p for g in base.k3_grad(x2))
    f2 = (base.k3_F(x2) // q2) % p
    coeffs = deep_coeffs_mod_p2(x2, p, pairs, deep_mask)
    return (
        all(g == 0 for g in grad),
        f2 == 0,
        coeff_signature(coeffs, p),
    )


def collect_p2_charts(p, samples_per_chart):
    q2 = p * p
    sc2 = base.squareclass_table(p, 2)
    pairs = {name: base.cross_pairs(left, right) for name, left, right in base.PARTITIONS}
    charts = {}
    ambiguous = modp_ambiguous_roots(p)
    t0 = time.time()

    for index, (x0, parts, sig, singular) in enumerate(ambiguous, 1):
        for lift in base.lift_vectors_mod_p2(x0, p):
            x2 = tuple(x0[i] + p * lift[i] for i in range(4))
            squares = [(v * v) % q2 for v in x2]
            for partition in parts:
                status, deep_mask, known = classify_partition_with_known(squares, q2, sc2, pairs[partition])
                if status != "deep":
                    continue
                key = broad_chart_key(p, partition, sig, deep_mask, x2, known, singular)
                if key not in charts:
                    charts[key] = {"count": 0, "samples": {}}
                charts[key]["count"] += 1
                if len(charts[key]["samples"]) < samples_per_chart:
                    ctrl = control_key(p, x2, pairs[partition], deep_mask)
                    if ctrl not in charts[key]["samples"]:
                        charts[key]["samples"][ctrl] = x2
        if index % 500 == 0:
            print(
                f"p={p}: collected {index}/{len(ambiguous)} ambiguous roots; "
                f"charts={len(charts)} elapsed={time.time() - t0:.1f}s",
                flush=True,
            )
    return charts


def p3_counts_for_sample(p, partition, x2, max_examples):
    q2 = p * p
    q3 = q2 * p
    sc3 = base.squareclass_table(p, 3)
    pairs = {name: base.cross_pairs(left, right) for name, left, right in base.PARTITIONS}
    counts = Counter()
    deep_masks = Counter()
    examples = {"killed": [], "deep": [], "resolved_ok": []}
    for lift in lift_vectors_mod_p3_from_p2(x2, p):
        x3 = tuple(x2[i] + q2 * lift[i] for i in range(4))
        squares = [(v * v) % q3 for v in x3]
        status, mask = base.classify_partition(squares, q3, sc3, pairs[partition])
        counts[status] += 1
        if status == "deep":
            deep_masks[mask] += 1
        if len(examples[status]) < max_examples:
            examples[status].append((x3, mask))
    return counts, deep_masks, examples


def outcome(counts):
    present = [status for status in STATUS_ORDER if counts[status]]
    return "+".join(present) if present else "no_p3_lift"


def analyze_prime(p, samples_per_chart, max_examples):
    t0 = time.time()
    charts = collect_p2_charts(p, samples_per_chart)
    summary = {
        "p": p,
        "chart_count": len(charts),
        "p2_deep_class_count": sum(item["count"] for item in charts.values()),
        "sample_count": sum(len(item["samples"]) for item in charts.values()),
        "chart_sample_outcomes_by_partition": defaultdict(Counter),
        "chart_union_outcomes_by_partition": defaultdict(Counter),
        "chart_union_outcomes_by_mask": defaultdict(Counter),
        "p2_weight_by_mask": defaultdict(int),
        "sample_p3_status_by_partition": defaultdict(Counter),
        "sample_p3_deep_masks_by_partition": defaultdict(Counter),
        "rows": [],
        "elapsed_seconds": 0.0,
    }

    for index, (key, item) in enumerate(charts.items(), 1):
        partition, sig, deep_mask, vp, rel, known, singular = key
        union = set()
        sample_rows = []
        aggregate_counts = Counter()
        aggregate_deep_masks = Counter()
        for ctrl, x2 in item["samples"].items():
            counts, deep_masks, examples = p3_counts_for_sample(p, partition, x2, max_examples)
            out = outcome(counts)
            union.add(out)
            aggregate_counts.update(counts)
            aggregate_deep_masks.update(deep_masks)
            summary["chart_sample_outcomes_by_partition"][partition][out] += 1
            summary["sample_p3_status_by_partition"][partition].update(counts)
            summary["sample_p3_deep_masks_by_partition"][partition].update(deep_masks)
            sample_rows.append({"control": ctrl, "rep": x2, "counts": counts, "deep_masks": deep_masks, "outcome": out, "examples": examples})

        union_key = ";".join(sorted(union)) if union else "no_samples"
        summary["chart_union_outcomes_by_partition"][partition][union_key] += 1
        summary["chart_union_outcomes_by_mask"][(partition, deep_mask)][union_key] += 1
        summary["p2_weight_by_mask"][(partition, deep_mask)] += item["count"]
        summary["rows"].append({
            "key": key,
            "p2_count": item["count"],
            "samples": sample_rows,
            "aggregate_counts": aggregate_counts,
            "aggregate_deep_masks": aggregate_deep_masks,
            "union_outcome": union_key,
        })

        if index % 50 == 0:
            print(
                f"p={p}: lifted samples for {index}/{len(charts)} charts "
                f"elapsed={time.time() - t0:.1f}s",
                flush=True,
            )

    summary["rows"].sort(key=lambda row: row["p2_count"], reverse=True)
    summary["elapsed_seconds"] = time.time() - t0
    return summary


def fmt_sig(sig):
    return "+".join(sig) if sig else "good"


def fmt_counts(counter):
    return " ".join(f"{status}={counter[status]}" for status in STATUS_ORDER)


def fmt_counter(counter):
    return " ".join(f"{key}={counter[key]}" for key in sorted(counter))


def write_report(path, summaries, top):
    with path.open("w") as out:
        out.write("A2244 sampled p^3 refinement of normalized p^2 deep boundary charts\n")
        out.write("all_zero_mod_p_root=excluded\n")
        out.write("scope=several p^3-control representatives per broad p^2 chart key\n")
        out.write("warning=diagnostic chart sampling, not exhaustive p^3 enumeration of all p^2 residues\n")
        out.write("status: deep means at least one relevant cross-difference is 0 mod p^3\n")

        for summary in summaries:
            p = summary["p"]
            out.write(f"\np={p}\n")
            out.write(f"elapsed_seconds={summary['elapsed_seconds']:.3f}\n")
            out.write(f"broad_p2_chart_keys={summary['chart_count']}\n")
            out.write(f"sampled_p2_representatives={summary['sample_count']}\n")
            out.write(f"p2_deep_root_partitions_represented={summary['p2_deep_class_count']}\n")

            out.write("\nsample outcomes by partition\n")
            for partition in sorted(summary["chart_sample_outcomes_by_partition"]):
                out.write(f"  {partition} {fmt_counter(summary['chart_sample_outcomes_by_partition'][partition])}\n")

            out.write("\nbroad chart union outcomes by partition\n")
            for partition in sorted(summary["chart_union_outcomes_by_partition"]):
                out.write(f"  {partition} {fmt_counter(summary['chart_union_outcomes_by_partition'][partition])}\n")

            out.write("\nsampled p^3 lift statuses by partition\n")
            for partition in sorted(summary["sample_p3_status_by_partition"]):
                out.write(f"  {partition} {fmt_counts(summary['sample_p3_status_by_partition'][partition])}\n")

            out.write("\nbroad chart union outcomes by p^2 deep mask\n")
            for (partition, mask), counter in sorted(summary["chart_union_outcomes_by_mask"].items(), key=lambda kv: (kv[0][0], kv[0][1])):
                out.write(
                    f"  {partition} {'+'.join(mask)} charts={sum(counter.values())} "
                    f"p2_weight={summary['p2_weight_by_mask'][(partition, mask)]} {fmt_counter(counter)}\n"
                )

            out.write("\ntop sampled p^3 deep masks by partition\n")
            for partition in sorted(summary["sample_p3_deep_masks_by_partition"]):
                out.write(f"  {partition}\n")
                for mask, count in summary["sample_p3_deep_masks_by_partition"][partition].most_common(top):
                    out.write(f"    {count} {'+'.join(mask)}\n")

            out.write("\ntop broad p^2 charts by represented class count\n")
            for row in summary["rows"][:top]:
                partition, sig, deep_mask, vp, rel, known, singular = row["key"]
                out.write(
                    f"  p2_count={row['p2_count']} partition={partition} boundary={fmt_sig(sig)} "
                    f"p2_deep={'+'.join(deep_mask)} vp={vp} rel={''.join(rel)} "
                    f"known={known} singular={singular} samples={len(row['samples'])} "
                    f"union={row['union_outcome']} aggregate_sample_counts=({fmt_counts(row['aggregate_counts'])})\n"
                )
                for sample in row["samples"][:3]:
                    out.write(
                        f"    sample control={sample['control']} rep={sample['rep']} "
                        f"outcome={sample['outcome']} counts=({fmt_counts(sample['counts'])})\n"
                    )

            out.write("\nexamples with resolved_ok sampled p^3 lifts\n")
            printed = 0
            for row in summary["rows"]:
                partition, sig, deep_mask, vp, rel, known, singular = row["key"]
                for sample in row["samples"]:
                    if not sample["counts"]["resolved_ok"]:
                        continue
                    example = sample["examples"]["resolved_ok"][0][0]
                    out.write(
                        f"  partition={partition} boundary={fmt_sig(sig)} p2_deep={'+'.join(deep_mask)} "
                        f"rep={sample['rep']} example_p3={example} counts=({fmt_counts(sample['counts'])})\n"
                    )
                    printed += 1
                    break
                if printed >= top:
                    break

            out.write("\nexamples with deep sampled p^3 lifts\n")
            printed = 0
            for row in summary["rows"]:
                partition, sig, deep_mask, vp, rel, known, singular = row["key"]
                for sample in row["samples"]:
                    if not sample["counts"]["deep"]:
                        continue
                    example, mask = sample["examples"]["deep"][0]
                    out.write(
                        f"  partition={partition} boundary={fmt_sig(sig)} p2_deep={'+'.join(deep_mask)} "
                        f"rep={sample['rep']} example_p3={example} p3_deep={'+'.join(mask)} "
                        f"counts=({fmt_counts(sample['counts'])})\n"
                    )
                    printed += 1
                    break
                if printed >= top:
                    break


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", nargs="+", type=int, default=[11, 23])
    parser.add_argument("--samples-per-chart", type=int, default=6)
    parser.add_argument("--max-examples", type=int, default=2)
    parser.add_argument("--top", type=int, default=16)
    parser.add_argument("--output", type=Path, default=Path("data/a2244_padic_deep_chart_p3_report.txt"))
    args = parser.parse_args()

    summaries = []
    for p in args.primes:
        print(f"analyzing p={p}", flush=True)
        summaries.append(analyze_prime(p, args.samples_per_chart, args.max_examples))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_report(args.output, summaries, args.top)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
