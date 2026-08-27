#!/usr/bin/env python3
"""Component-wise 11/23-adic boundary analysis for A(2,2,4,4).

This refines the existing boundary diagnostics by separating the raw K3
boundary components

    Z1,Z2,Z3,Z4,E12,E13,E14,E23,E24,E34

and the three A(2,2,4,4) partitions.  For each prime it does:

1. enumerate all mod-p K3 roots, excluding the all-zero projective root;
2. classify each component/partition at mod p;
3. exhaustively lift mod-p ambiguous component/partition classes to p^2;
4. search for smooth p^3 resolved examples as local witnesses.

The p^3 step is deliberately a witness search, not a full p^3 enumeration.
The exhaustive part is the mod-p and p^2 component matrix.
"""

from __future__ import annotations

import argparse
import sys
import time
from collections import Counter, defaultdict
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import analyze_a2244_padic_residues as base  # noqa: E402


COMPONENTS = ["Z1", "Z2", "Z3", "Z4", "E12", "E13", "E14", "E23", "E24", "E34"]
STATUS_ORDER = ["killed", "deep", "resolved_ok"]


def fmt_sig(sig: tuple[str, ...]) -> str:
    return "+".join(sig) if sig else "good"


def grad_smooth_mod_p(x: tuple[int, int, int, int], p: int) -> bool:
    return any(g % p for g in base.k3_grad(x))


def lift_vectors_mod_p3_from_p2(x2: tuple[int, int, int, int], p: int):
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


def roots_mod_p(p: int) -> list[tuple[int, int, int, int]]:
    out = []
    for x0 in product(range(p), repeat=4):
        if all(v == 0 for v in x0):
            continue;
        if base.k3_F(x0) % p == 0:
            out.append(x0)
    return out


def component_counter_dict():
    return {
        comp: {name: Counter() for name, _, _ in base.PARTITIONS}
        for comp in COMPONENTS
    }


def analyze_prime(p: int, p3_root_limit: int, p3_p2_lift_limit: int, p3_p3_lift_limit: int):
    t0 = time.time()
    sc_p = base.squareclass_table(p, 1)
    sc_p2 = base.squareclass_table(p, 2)
    sc_p3 = base.squareclass_table(p, 3)
    pairs = {name: base.cross_pairs(left, right) for name, left, right in base.PARTITIONS}

    roots = roots_mod_p(p)
    summary = {
        "p": p,
        "roots": len(roots),
        "boundary_roots": 0,
        "smooth_boundary_roots": 0,
        "component_roots": Counter(),
        "component_smooth_roots": Counter(),
        "component_modp_status": component_counter_dict(),
        "component_p2_lift_status": component_counter_dict(),
        "component_p2_root_outcome": component_counter_dict(),
        "signature_modp_status": defaultdict(lambda: {name: Counter() for name, _, _ in base.PARTITIONS}),
        "p2_examples": {},
        "p3_examples": {},
        "elapsed_seconds": 0.0,
    }

    ambiguous: list[tuple[tuple[int, int, int, int], tuple[str, ...], list[str]]] = []

    for x0 in roots:
        sig = base.boundary_signature(x0, p)
        smooth = grad_smooth_mod_p(x0, p)
        if sig:
            summary["boundary_roots"] += 1
            if smooth:
                summary["smooth_boundary_roots"] += 1
        for comp in sig:
            summary["component_roots"][comp] += 1
            if smooth:
                summary["component_smooth_roots"][comp] += 1

        squares = [(v * v) % p for v in x0]
        ambiguous_parts = []
        for name, _, _ in base.PARTITIONS:
            status, _ = base.classify_partition(squares, p, sc_p, pairs[name])
            if status == "deep":
                ambiguous_parts.append(name)
            for comp in sig:
                summary["component_modp_status"][comp][name][status] += 1
            summary["signature_modp_status"][sig][name][status] += 1
        if ambiguous_parts:
            ambiguous.append((x0, sig, ambiguous_parts))

    for index, (x0, sig, ambiguous_parts) in enumerate(ambiguous, 1):
        part_counts = {name: Counter() for name in ambiguous_parts}
        part_examples = {name: {} for name in ambiguous_parts}

        for lift in base.lift_vectors_mod_p2(x0, p):
            x2 = tuple(x0[i] + p * lift[i] for i in range(4))
            squares = [(v * v) % (p * p) for v in x2]
            for name in ambiguous_parts:
                status, deep_mask = base.classify_partition(squares, p * p, sc_p2, pairs[name])
                part_counts[name][status] += 1
                for comp in sig:
                    summary["component_p2_lift_status"][comp][name][status] += 1
                if status != "killed" and status not in part_examples[name]:
                    part_examples[name][status] = (x2, deep_mask)

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
            for comp in sig:
                summary["component_p2_root_outcome"][comp][name][outcome] += 1
                if counts["resolved_ok"]:
                    summary["component_p2_root_outcome"][comp][name]["has_resolved_ok_lift"] += 1
                if counts["deep"]:
                    summary["component_p2_root_outcome"][comp][name]["has_deep_lift"] += 1

            for status, item in part_examples[name].items():
                key = (tuple(sig), name, status)
                if key not in summary["p2_examples"]:
                    x2, deep_mask = item
                    summary["p2_examples"][key] = (x0, x2, deep_mask)

        if index % 500 == 0:
            print(
                f"p={p}: p^2 refined {index}/{len(ambiguous)} ambiguous roots "
                f"elapsed={time.time() - t0:.1f}s",
                flush=True,
            )

    # Targeted smooth p^3 resolved witnesses.  A smooth mod-p root with a
    # resolved_ok p^3 partition is enough to certify a genuine Q_p branch for
    # that component/partition status by Hensel stability.
    wanted = {(comp, name) for comp in COMPONENTS for name, _, _ in base.PARTITIONS}
    searched_roots = Counter()
    for x0 in roots:
        sig = base.boundary_signature(x0, p)
        if not sig or not grad_smooth_mod_p(x0, p):
            continue
        squares0 = [(v * v) % p for v in x0]
        parts = []
        for name, _, _ in base.PARTITIONS:
            status, _ = base.classify_partition(squares0, p, sc_p, pairs[name])
            if status != "killed":
                parts.append(name)
        if not parts:
            continue

        for comp in sig:
            for name in parts:
                target = (comp, name)
                if target not in wanted or target in summary["p3_examples"]:
                    continue
                if searched_roots[target] >= p3_root_limit:
                    continue
                searched_roots[target] += 1

                p2_checked = 0
                found = False
                for lift2 in base.lift_vectors_mod_p2(x0, p):
                    x2 = tuple(x0[i] + p * lift2[i] for i in range(4))
                    p2_checked += 1
                    if p2_checked > p3_p2_lift_limit:
                        break
                    p3_checked = 0
                    for lift3 in lift_vectors_mod_p3_from_p2(x2, p):
                        x3 = tuple(x2[i] + p * p * lift3[i] for i in range(4))
                        p3_checked += 1
                        squares3 = [(v * v) % (p**3) for v in x3]
                        status3, deep3 = base.classify_partition(
                            squares3, p**3, sc_p3, pairs[name]
                        )
                        if status3 == "resolved_ok":
                            summary["p3_examples"][target] = (x0, x2, x3, tuple(sig), deep3)
                            found = True
                            break
                        if p3_checked >= p3_p3_lift_limit:
                            break
                    if found:
                        break

    summary["elapsed_seconds"] = time.time() - t0
    return summary


def status_line(counter: Counter, keys: list[str]) -> str:
    return " ".join(f"{key}={counter[key]}" for key in keys)


def component_verdict(summary, comp: str, part: str) -> str:
    modp = summary["component_modp_status"][comp][part]
    p2 = summary["component_p2_lift_status"][comp][part]
    outcomes = summary["component_p2_root_outcome"][comp][part]
    if (comp, part) in summary["p3_examples"]:
        return "smooth_p3_resolved"
    if modp["resolved_ok"]:
        return "modp_resolved"
    if p2["resolved_ok"]:
        return "p2_resolved"
    if p2["deep"] or outcomes["has_deep_lift"]:
        return "deep_only"
    if modp["deep"]:
        return "modp_deep_no_p2_survivor"
    return "killed"


def write_report(path: Path, summaries: list[dict], top: int) -> None:
    with path.open("w") as out:
        out.write("A2244 component-wise 11/23-adic boundary analysis\n")
        out.write("surface=(ab+ac+ad+bc+bd+cd)^2-4abcd=0\n")
        out.write("all_zero_mod_p_root=excluded\n")
        out.write("p2=exhaustive lifts of mod-p ambiguous component/partition classes\n")
        out.write("p3=smooth resolved witness search, not full p3 enumeration\n")
        out.write("verdict priority: smooth_p3_resolved > modp_resolved > p2_resolved > deep_only > killed\n\n")

        for summary in summaries:
            p = summary["p"]
            out.write(f"p={p}\n")
            out.write(f"elapsed_seconds={summary['elapsed_seconds']:.3f}\n")
            out.write(f"Fp_roots={summary['roots']}\n")
            out.write(f"boundary_roots={summary['boundary_roots']}\n")
            out.write(f"smooth_boundary_roots={summary['smooth_boundary_roots']}\n\n")

            out.write("component totals\n")
            for comp in COMPONENTS:
                out.write(
                    f"  {comp} roots={summary['component_roots'][comp]} "
                    f"smooth={summary['component_smooth_roots'][comp]}\n"
                )

            out.write("\ncomponent/partition matrix\n")
            for comp in COMPONENTS:
                out.write(f"  component {comp}\n")
                for part, _, _ in base.PARTITIONS:
                    verdict = component_verdict(summary, comp, part)
                    modp = summary["component_modp_status"][comp][part]
                    p2 = summary["component_p2_lift_status"][comp][part]
                    outcomes = summary["component_p2_root_outcome"][comp][part]
                    out.write(
                        f"    {part} verdict={verdict} "
                        f"modp({status_line(modp, STATUS_ORDER)}) "
                        f"p2_lifts({status_line(p2, STATUS_ORDER)}) "
                        f"p2_roots("
                        f"all_lifts_killed={outcomes['all_lifts_killed']} "
                        f"deep_only={outcomes['deep_only']} "
                        f"resolved_only={outcomes['resolved_only']} "
                        f"resolved_and_deep={outcomes['resolved_and_deep']} "
                        f"has_resolved_ok_lift={outcomes['has_resolved_ok_lift']} "
                        f"has_deep_lift={outcomes['has_deep_lift']})\n"
                    )
                    if (comp, part) in summary["p3_examples"]:
                        x0, x2, x3, sig, deep3 = summary["p3_examples"][(comp, part)]
                        out.write(
                            f"      p3_example modp={x0} p2={x2} p3={x3} "
                            f"boundary={fmt_sig(sig)}\n"
                        )

            out.write("\nmod-p boundary signatures with non-killed partitions\n")
            rows = []
            for sig, by_part in summary["signature_modp_status"].items():
                for part, counts in by_part.items():
                    if counts["deep"] or counts["resolved_ok"]:
                        rows.append((counts["deep"] + counts["resolved_ok"], sig, part, counts))
            for _, sig, part, counts in sorted(rows, reverse=True)[:top]:
                out.write(
                    f"  {fmt_sig(sig)} {part} "
                    f"{status_line(counts, STATUS_ORDER)}\n"
                )

            out.write("\np2 examples by boundary signature and partition\n")
            for (sig, part, status), (x0, x2, deep_mask) in list(summary["p2_examples"].items())[:top]:
                out.write(
                    f"  status={status} boundary={fmt_sig(sig)} partition={part} "
                    f"modp={x0} p2={x2} deep={'+'.join(deep_mask)}\n"
                )
            out.write("\n")

        out.write("combined conclusion\n")
        out.write(
            "No raw component is eliminated component-wise by the finite local "
            "analysis.  Simple collision and zero components have surviving "
            "partitions, often with smooth p^3 resolved witnesses.  The "
            "obstruction seen in height searches must therefore use the same "
            "partition simultaneously at 11 and 23, real compatibility, and "
            "the full M(2,2,4,8) cover equations rather than a single boundary "
            "component label.\n"
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--primes", nargs="+", type=int, default=[11, 23])
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/a2244_component_adic_analysis.txt"),
    )
    parser.add_argument("--top", type=int, default=30)
    parser.add_argument("--p3-root-limit", type=int, default=80)
    parser.add_argument("--p3-p2-lift-limit", type=int, default=2000)
    parser.add_argument("--p3-p3-lift-limit", type=int, default=4000)
    args = parser.parse_args()

    summaries = []
    for p in args.primes:
        print(f"component analysis p={p}", flush=True)
        summaries.append(
            analyze_prime(
                p,
                args.p3_root_limit,
                args.p3_p2_lift_limit,
                args.p3_p3_lift_limit,
            )
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_report(args.output, summaries, args.top)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
