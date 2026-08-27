#!/usr/bin/env python3
"""Candidate-level p-adic signature sieve for A(2,2,4,4).

This is the finite candidate sieve suggested by the p^2/p^3 boundary chart
analysis.  It reads integer K3 surface tuples, computes the oriented
A(2,2,4,4) cross-difference squareclass status for each 2+2 partition at
11^k and 23^k, and keeps only tuples with a common finite-local partition.

The test is intentionally signed: over Q_p the squareclass of -u matters.
This is important at p=11 and p=23, where -1 is nonsquare.
"""

from __future__ import annotations

import argparse
from collections import Counter
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


def parse_tuple(line: str) -> tuple[int, int, int, int] | None:
    line = line.strip()
    if not (line.startswith("[") and line.endswith("]")):
        return None
    vals = tuple(int(x) for x in line[1:-1].split(","))
    return vals if len(vals) == 4 else None


def read_tuples(path: Path) -> list[tuple[int, int, int, int]]:
    return [
        tup
        for line in path.read_text().splitlines()
        for tup in [parse_tuple(line)]
        if tup is not None
    ]


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
        raise ValueError("Legendre symbol requested for zero")
    return 1 if pow(u, (p - 1) // 2, p) == 1 else -1


def vp_unit_signed(n: int, p: int) -> tuple[int | None, int]:
    if n == 0:
        return None, 0
    e = 0
    while n % p == 0:
        e += 1
        n //= p
    return e, n % p


def squareclass_qp_signed(n: int, p: int) -> tuple[int, int] | str:
    if n == 0:
        return "0"
    e, unit = vp_unit_signed(n, p)
    assert e is not None
    return (e % 2, legendre(unit, p))


def squareclass_mod_pk(n: int, p: int, k: int) -> tuple[int, int] | None:
    q = p**k
    r = n % q
    if r == 0:
        return None
    e = 0
    while r % p == 0:
        e += 1
        r //= p
    return (e % 2, legendre(r, p))


def cross_pairs(left: tuple[int, int], right: tuple[int, int]) -> list[tuple[int, int]]:
    return [(i, j) for i in left for j in right]


def pair_label(i: int, j: int) -> str:
    return PAIR_LABELS[(i, j) if i < j else (j, i)]


def boundary_signature(t: tuple[int, int, int, int], p: int) -> tuple[str, ...]:
    labels: list[str] = []
    residues = [x % p for x in t]
    squares = [(x * x) % p for x in residues]
    for i, x in enumerate(residues):
        if x == 0:
            labels.append(f"Z{i + 1}")
    for i in range(4):
        for j in range(i + 1, 4):
            if squares[i] == squares[j]:
                labels.append(PAIR_LABELS[(i, j)])
    return tuple(labels)


def partition_status_mod_pk(
    t: tuple[int, int, int, int],
    p: int,
    k: int,
    left: tuple[int, int],
    right: tuple[int, int],
) -> tuple[str, tuple[tuple[int, int], ...], tuple[str, ...]]:
    q = p**k
    squares = [(x * x) % q for x in t]
    known: set[tuple[int, int]] = set()
    deep: list[str] = []
    for i, j in cross_pairs(left, right):
        sc = squareclass_mod_pk(squares[i] - squares[j], p, k)
        if sc is None:
            deep.append(pair_label(i, j))
        else:
            known.add(sc)
    if len(known) > 1:
        return "killed", tuple(sorted(known)), ()
    if deep:
        return "deep", tuple(sorted(known)), tuple(deep)
    return "resolved_ok", tuple(sorted(known)), ()


def partition_status_exact_qp(
    t: tuple[int, int, int, int],
    p: int,
    left: tuple[int, int],
    right: tuple[int, int],
) -> tuple[bool, tuple[tuple[int, int] | str, ...]]:
    squares = [x * x for x in t]
    classes = tuple(
        squareclass_qp_signed(squares[i] - squares[j], p)
        for i, j in cross_pairs(left, right)
    )
    return len(set(classes)) == 1, classes


def possible_partitions_mod_pk(
    t: tuple[int, int, int, int], p: int, k: int
) -> tuple[str, ...]:
    out = []
    for name, left, right in PARTITIONS:
        status, _, _ = partition_status_mod_pk(t, p, k, left, right)
        if status != "killed":
            out.append(name)
    return tuple(out)


def exact_qp_partitions(t: tuple[int, int, int, int], p: int) -> tuple[str, ...]:
    out = []
    for name, left, right in PARTITIONS:
        ok, _ = partition_status_exact_qp(t, p, left, right)
        if ok:
            out.append(name)
    return tuple(out)


def real_partitions_for_sorted_positive(t: tuple[int, int, int, int]) -> tuple[str, ...]:
    if tuple(t) != tuple(sorted(t)) or any(x <= 0 for x in t):
        return real_partitions_by_sign(t)
    return ("12|34",)


def real_partitions_by_sign(t: tuple[int, int, int, int]) -> tuple[str, ...]:
    squares = [x * x for x in t]
    out = []
    for name, left, right in PARTITIONS:
        products = []
        for i in left:
            products.append((squares[i] - squares[right[0]]) * (squares[i] - squares[right[1]]))
        for j in right:
            products.append((squares[left[0]] - squares[j]) * (squares[left[1]] - squares[j]))
        if all(x >= 0 for x in products):
            out.append(name)
    return tuple(out)


def k3_p3_lift_control(t: tuple[int, int, int, int], p: int) -> str:
    q2 = p * p
    x2 = tuple(x % q2 for x in t)
    f2 = k3_F(x2)
    grad = tuple(g % p for g in k3_grad(x2))
    if any(grad):
        return "smooth"
    return "singular_lift" if f2 % (p**3) == 0 else "singular_no_lift"


def deep_coeff_signature(
    t: tuple[int, int, int, int],
    p: int,
    left: tuple[int, int],
    right: tuple[int, int],
    deep_mask: tuple[str, ...],
) -> tuple[tuple[str, int | str], ...]:
    q2 = p * p
    q3 = q2 * p
    residues = [x % q3 for x in t]
    deep = set(deep_mask)
    out = []
    for i, j in cross_pairs(left, right):
        label = pair_label(i, j)
        if label not in deep:
            continue
        diff = (residues[i] * residues[i] - residues[j] * residues[j]) % q3
        coeff = (diff // q2) % p
        if coeff == 0:
            out.append((label, 0))
        else:
            out.append((label, legendre(coeff, p)))
    return tuple(out)


def signature_for_partition(
    t: tuple[int, int, int, int],
    p: int,
    k: int,
    name: str,
    left: tuple[int, int],
    right: tuple[int, int],
) -> tuple:
    status, known, deep = partition_status_mod_pk(t, p, k, left, right)
    _, _, deep2 = partition_status_mod_pk(t, p, 2, left, right)
    return (
        name,
        boundary_signature(t, p),
        k3_p3_lift_control(t, p),
        status,
        known,
        deep,
        deep2,
        deep_coeff_signature(t, p, left, right, deep2) if deep2 else (),
    )


def summarize(tuples: list[tuple[int, int, int, int]], primes: tuple[int, int], max_depth: int):
    summary = {
        "by_depth": {},
        "depth_distribution": {k: Counter() for k in range(1, max_depth + 1)},
        "p3_signature_counts": {p: Counter() for p in primes},
        "p2_finite_common": [],
        "p3_finite_common": [],
        "exact_signed_finite_common": [],
    }
    real_parts_by_tuple = {t: set(real_partitions_for_sorted_positive(t)) for t in tuples}

    for k in range(1, max_depth + 1):
        depth_rows = []
        for t in tuples:
            possible = {p: set(possible_partitions_mod_pk(t, p, k)) for p in primes}
            common = set.intersection(*(possible[p] for p in primes))
            real_common = common & real_parts_by_tuple[t]
            key = (
                tuple(sorted(possible[primes[0]])),
                tuple(sorted(possible[primes[1]])),
                tuple(sorted(common)),
                tuple(sorted(real_common)),
            )
            summary["depth_distribution"][k][key] += 1
            if common:
                depth_rows.append((t, tuple(sorted(common)), tuple(sorted(real_common)), possible))
        summary["by_depth"][k] = depth_rows

    summary["p2_finite_common"] = summary["by_depth"].get(2, [])
    summary["p3_finite_common"] = summary["by_depth"].get(3, [])

    for t in tuples:
        exact = {p: set(exact_qp_partitions(t, p)) for p in primes}
        common = set.intersection(*(exact[p] for p in primes))
        real_common = common & real_parts_by_tuple[t]
        if common:
            summary["exact_signed_finite_common"].append(
                (t, tuple(sorted(common)), tuple(sorted(real_common)), exact)
            )
        for p in primes:
            for name, left, right in PARTITIONS:
                sig = signature_for_partition(t, p, 3, name, left, right)
                summary["p3_signature_counts"][p][sig] += 1
    return summary


def write_tuple_file(path: Path, rows) -> None:
    with path.open("w") as out:
        for t, *_ in rows:
            out.write("[" + ",".join(str(x) for x in t) + "]\n")


def write_report(path: Path, input_path: Path, tuples, primes, summary, top: int) -> None:
    with path.open("w") as out:
        out.write("A2244 p-adic signature sieve for K3 surface tuples\n")
        out.write(f"input={input_path}\n")
        out.write(f"tuples={len(tuples)}\n")
        out.write(f"primes={','.join(str(p) for p in primes)}\n")
        out.write("squareclasses=signed\n")
        out.write("possible_statuses=resolved_ok_or_deep\n")
        out.write("k3_lift_control=positive-coordinate diagnostic only; do not use as a killer without a signed surface representative\n")

        out.write("\nDepth summary\n")
        for k, rows in summary["by_depth"].items():
            real_count = sum(1 for _, _, real_common, _ in rows if real_common)
            out.write(
                f"k={k} finite_common={len(rows)} finite_common_and_real={real_count}\n"
            )
            for key, count in summary["depth_distribution"][k].most_common(10):
                out.write(f"  {count} {key}\n")

        out.write("\nExact signed Qp summary\n")
        exact_rows = summary["exact_signed_finite_common"]
        real_exact = sum(1 for _, _, real_common, _ in exact_rows if real_common)
        out.write(
            f"finite_common={len(exact_rows)} finite_common_and_real={real_exact}\n"
        )
        if exact_rows:
            for t, common, real_common, exact in exact_rows[:top]:
                out.write(
                    f"  tuple={t} common={common} real_common={real_common} "
                    f"p{primes[0]}={tuple(sorted(exact[primes[0]]))} "
                    f"p{primes[1]}={tuple(sorted(exact[primes[1]]))}\n"
                )

        out.write("\np^2 finite-common tuples and p^3 disposition\n")
        p3_survivor_tuples = {row[0] for row in summary["p3_finite_common"]}
        for t, common, real_common, possible in summary["p2_finite_common"]:
            p3_parts = {
                p: set(possible_partitions_mod_pk(t, p, 3)) for p in primes
            }
            p3_common = tuple(sorted(set.intersection(*(p3_parts[p] for p in primes))))
            out.write(
                f"tuple={t} p2_common={common} p2_real={real_common} "
                f"p3_common={p3_common} disposition="
                f"{'survives_p3' if t in p3_survivor_tuples else 'killed_by_p3'}\n"
            )
            for p in primes:
                out.write(f"  p={p}\n")
                for name, left, right in PARTITIONS:
                    sig = signature_for_partition(t, p, 3, name, left, right)
                    out.write(f"    {sig}\n")

        out.write("\np^3 finite-common tuples\n")
        if not summary["p3_finite_common"]:
            out.write("none\n")
        for t, common, real_common, possible in summary["p3_finite_common"][:top]:
            out.write(
                f"tuple={t} common={common} real_common={real_common} "
                f"p{primes[0]}={tuple(sorted(possible[primes[0]]))} "
                f"p{primes[1]}={tuple(sorted(possible[primes[1]]))}\n"
            )

        out.write("\nTop p^3 signatures by prime\n")
        for p in primes:
            out.write(f"p={p}\n")
            for sig, count in summary["p3_signature_counts"][p].most_common(top):
                out.write(f"  {count} {sig}\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input",
        nargs="?",
        type=Path,
        default=Path("data/surface_tuples_B10000_strata_non_d_zero_11_23.txt"),
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/a2244_padic_signature_sieve_B10000.txt"),
    )
    parser.add_argument(
        "--survivors",
        type=Path,
        default=Path("data/a2244_padic_signature_sieve_B10000_p3_survivors.txt"),
    )
    parser.add_argument("--primes", nargs=2, type=int, default=[11, 23])
    parser.add_argument("--max-depth", type=int, default=3)
    parser.add_argument("--top", type=int, default=20)
    args = parser.parse_args()

    tuples = read_tuples(args.input)
    primes = tuple(args.primes)
    summary = summarize(tuples, primes, args.max_depth)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    write_report(args.output, args.input, tuples, primes, summary, args.top)
    final_rows = summary["by_depth"].get(args.max_depth, [])
    write_tuple_file(args.survivors, final_rows)

    print(f"wrote {args.output}")
    print(f"wrote {args.survivors}")
    print(f"p^{args.max_depth} finite-common survivors: {len(final_rows)}")


if __name__ == "__main__":
    main()
