#!/usr/bin/env python3
from collections import Counter
from itertools import permutations
from pathlib import Path
import sys

PARTITIONS = [
    ((0, 1), (2, 3), "12|34"),
    ((0, 2), (1, 3), "13|24"),
    ((0, 3), (1, 2), "14|23"),
]
PRIMES = [11, 23]
AUX_PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def parse_tuple(line):
    line = line.strip()
    if not (line.startswith("[") and line.endswith("]")):
        return None
    vals = tuple(int(x) for x in line[1:-1].split(","))
    return vals if len(vals) == 4 else None


def read_tuples(path):
    return [t for line in Path(path).read_text().splitlines() for t in [parse_tuple(line)] if t]


def vp_unit(n, p):
    if n == 0:
        return None, 0
    n = int(n)
    e = 0
    while n % p == 0:
        e += 1
        n //= p
    return e, n % p


def is_qp_square(n, p):
    if n == 0:
        return True
    e, u = vp_unit(n, p)
    return e % 2 == 0 and pow(u, (p - 1) // 2, p) == 1


def squareclass_qp(n, p):
    if n == 0:
        return "0"
    e, u = vp_unit(n, p)
    leg = 1 if pow(u, (p - 1) // 2, p) == 1 else -1
    return (e % 2, leg)


def products_for_ordered(t):
    a, b, c, d = t
    A, B, C, D = [x * x for x in (a, b, c, d)]
    return [
        (A - C) * (A - D),
        (B - C) * (B - D),
        (A - C) * (B - C),
        (A - D) * (B - D),
    ]


def local_ok_ordered(t, p):
    return all(is_qp_square(v, p) for v in products_for_ordered(t))


def local_ok_same_ordering(t, primes):
    for perm in permutations(range(4)):
        ot = tuple(t[i] for i in perm)
        if all(local_ok_ordered(ot, p) for p in primes):
            return True
    return False


def real_ok_ordered(t):
    return all(v >= 0 for v in products_for_ordered(t))


def local_and_real_ok_same_ordering(t, primes):
    for perm in permutations(range(4)):
        ot = tuple(t[i] for i in perm)
        if real_ok_ordered(ot) and all(local_ok_ordered(ot, p) for p in primes):
            return True
    return False


def partition_ok(t, p, part):
    (i, j), (k, l), name = part
    X = [x * x for x in t]
    cross = [X[i] - X[k], X[i] - X[l], X[j] - X[k], X[j] - X[l]]
    classes = tuple(squareclass_qp(v, p) for v in cross)
    return len(set(classes)) == 1, classes


def ok_partitions(t, p):
    return tuple(name for part in PARTITIONS if partition_ok(t, p, part)[0] for name in [part[2]])


def boundary_components(t, p):
    labels = []
    for i, x in enumerate(t):
        if x % p == 0:
            labels.append(f"Z{i+1}")
    sq = [(x * x) % p for x in t]
    for i in range(4):
        for j in range(i + 1, 4):
            if sq[i] == sq[j]:
                labels.append(f"E{i+1}{j+1}")
    return tuple(labels)


def main():
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/surface_tuples_B10000_strata_non_d_zero_11_23.txt")
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/a2244_boundary_strata_diagnostics_B10000.txt")
    tuples = read_tuples(input_path)

    qlocal_summary = Counter()
    partition_distribution = Counter()
    signature_counts = {p: Counter() for p in PRIMES}
    signature_ok = {p: Counter() for p in PRIMES}
    component_counts = {p: Counter() for p in PRIMES}
    component_ok = {p: Counter() for p in PRIMES}
    common_survivors = []

    for t in tuples:
        ok_by_prime = {p: bool(ok_partitions(t, p)) for p in PRIMES}
        common_parts = tuple(sorted(set(ok_partitions(t, 11)) & set(ok_partitions(t, 23))))
        same_ordering = local_ok_same_ordering(t, PRIMES)
        same_ordering_real = local_and_real_ok_same_ordering(t, PRIMES)

        if same_ordering:
            qlocal_summary["same_ordering_Q11_Q23_local"] += 1
            common_survivors.append(t)
        else:
            qlocal_summary["no_same_ordering_Q11_Q23_local"] += 1
        if same_ordering_real:
            qlocal_summary["same_ordering_Q11_Q23_and_real"] += 1
        else:
            qlocal_summary["no_same_ordering_Q11_Q23_and_real"] += 1
        for p in PRIMES:
            if ok_by_prime[p]:
                qlocal_summary[f"Q{p}_local_possible"] += 1
            else:
                qlocal_summary[f"Q{p}_local_obstructed"] += 1

            comps = boundary_components(t, p)
            signature_counts[p][comps] += 1
            if ok_by_prime[p]:
                signature_ok[p][comps] += 1
            for comp in comps:
                component_counts[p][comp] += 1
                if ok_by_prime[p]:
                    component_ok[p][comp] += 1

        partition_distribution[(ok_partitions(t, 11), ok_partitions(t, 23), common_parts)] += 1

    with output_path.open("w") as out:
        out.write("A2244 boundary-stratum p-adic diagnostics\n")
        out.write(f"input={input_path}\n")
        out.write(f"tuples={len(tuples)}\n")
        out.write("\nQp local summary\n")
        for key in sorted(qlocal_summary):
            out.write(f"{key} {qlocal_summary[key]}\n")

        out.write("\nSurviving partition-set distribution\n")
        for key, count in partition_distribution.most_common():
            out.write(f"{count} {key}\n")

        out.write("\nPer-prime boundary signature totals and Qp-local counts\n")
        for p in PRIMES:
            out.write(f"p={p}\n")
            for sig, count in signature_counts[p].most_common():
                out.write(f"  {count} {signature_ok[p][sig]} {'+'.join(sig)}\n")

        out.write("\nPer-prime component marginal totals and Qp-local counts\n")
        for p in PRIMES:
            out.write(f"p={p}\n")
            for comp, count in component_counts[p].most_common():
                out.write(f"  {comp} {count} {component_ok[p][comp]}\n")

        out.write("\nSame-ordering Q11/Q23 local survivors\n")
        if not common_survivors:
            out.write("none\n")
        for t in common_survivors:
            out.write(f"tuple={t}\n")
            out.write(f"  boundary_11={boundary_components(t, 11)}\n")
            out.write(f"  boundary_23={boundary_components(t, 23)}\n")
            out.write(f"  ok_partitions_11={ok_partitions(t, 11)}\n")
            out.write(f"  ok_partitions_23={ok_partitions(t, 23)}\n")
            for perm in permutations(range(4)):
                ot = tuple(t[i] for i in perm)
                if all(local_ok_ordered(ot, p) for p in PRIMES):
                    bad = [p for p in AUX_PRIMES if not local_ok_ordered(ot, p)]
                    out.write(f"  perm={perm} ordered={ot} real_ok={real_ok_ordered(ot)} finite_bad_primes={bad}\n")
                    out.write(f"    products={products_for_ordered(ot)}\n")

    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
