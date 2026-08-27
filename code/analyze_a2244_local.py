#!/usr/bin/env python3
from itertools import combinations
from pathlib import Path
import sys

PARTITIONS = [((0, 1), (2, 3)), ((0, 2), (1, 3)), ((0, 3), (1, 2))]
PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47]


def parse_tuple_line(line):
    line = line.strip()
    if not line.startswith("[") or not line.endswith("]"):
        return None
    return [int(x) for x in line[1:-1].split(",")]


def read_tuples(path):
    out = []
    for line in Path(path).read_text().splitlines():
        tup = parse_tuple_line(line)
        if tup is not None and len(tup) == 4:
            out.append(tup)
    return out


def is_square_mod(x, p):
    x %= p
    if x == 0:
        return True
    return pow(x, (p - 1) // 2, p) == 1


def chi(x, p):
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def good_reduction(tup, p):
    if any(x % p == 0 for x in tup):
        return False
    squares = [(x * x) % p for x in tup]
    return len(set(squares)) == 4


def a2244_possible_original(tup, p):
    xs = [(x * x) % p for x in tup]
    for perm in ((0,1,2,3), (0,1,3,2), (0,2,1,3), (0,2,3,1), (0,3,1,2), (0,3,2,1),
                 (1,0,2,3), (1,0,3,2), (1,2,0,3), (1,2,3,0), (1,3,0,2), (1,3,2,0),
                 (2,0,1,3), (2,0,3,1), (2,1,0,3), (2,1,3,0), (2,3,0,1), (2,3,1,0),
                 (3,0,1,2), (3,0,2,1), (3,1,0,2), (3,1,2,0), (3,2,0,1), (3,2,1,0)):
        A, B, C, D = [xs[i] for i in perm]
        conditions = [
            (A - C) * (A - D),
            (B - C) * (B - D),
            (A - C) * (B - C),
            (A - D) * (B - D),
        ]
        if all(is_square_mod(v, p) for v in conditions):
            return True
    return False


def partition_patterns(tup, p):
    xs = [(x * x) % p for x in tup]
    rows = []
    for left, right in PARTITIONS:
        i, j = left
        k, l = right
        pattern = [
            chi(xs[i] - xs[k], p),
            chi(xs[i] - xs[l], p),
            chi(xs[j] - xs[k], p),
            chi(xs[j] - xs[l], p),
        ]
        rows.append((left, right, pattern))
    return rows


def a2244_possible_good_criterion(tup, p):
    for _, _, pattern in partition_patterns(tup, p):
        if 0 not in pattern and len(set(pattern)) == 1:
            return True
    return False


def surface_eq_mod(vals, p):
    a, b, c, d = [x % p for x in vals]
    s2 = (a*b + a*c + a*d + b*c + b*d + c*d) % p
    s4 = (a*b*c*d) % p
    return (s2*s2 - 4*s4) % p == 0


def finite_field_surface_counts(primes):
    rows = []
    for p in primes:
        signed_total = 0
        signed_good = 0
        signed_good_possible = 0
        class_all = set()
        class_good = set()
        class_good_possible = set()
        vals = range(1, p)
        for a in vals:
            for b in vals:
                for c in vals:
                    # Solve for d by direct loop. The largest p here is small.
                    for d in vals:
                        if not surface_eq_mod((a, b, c, d), p):
                            continue
                        signed_total += 1
                        xs = tuple(sorted([(a*a) % p, (b*b) % p, (c*c) % p, (d*d) % p]))
                        class_all.add(xs)
                        tup = [a, b, c, d]
                        if good_reduction(tup, p):
                            signed_good += 1
                            class_good.add(xs)
                            if a2244_possible_good_criterion(tup, p):
                                signed_good_possible += 1
                                class_good_possible.add(xs)
        rows.append({
            "p": p,
            "signed_total": signed_total,
            "signed_good": signed_good,
            "signed_good_possible": signed_good_possible,
            "classes_total": len(class_all),
            "classes_good": len(class_good),
            "classes_good_possible": len(class_good_possible),
        })
    return rows


def minimal_cover(obstructs, universe):
    primes = list(obstructs)
    for r in range(1, len(primes) + 1):
        solutions = []
        for combo in combinations(primes, r):
            covered = set()
            for p in combo:
                covered |= obstructs[p]
            if covered >= universe:
                solutions.append(combo)
        if solutions:
            return solutions
    return []


def main():
    input_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("data/surface_tuples_B2000.txt")
    output_path = Path(sys.argv[2]) if len(sys.argv) > 2 else Path("data/a2244_local_obstruction_analysis_B2000.txt")
    tuples = read_tuples(input_path)
    universe = set(range(len(tuples)))

    obstructs = {p: set() for p in PRIMES}
    possible = {p: set() for p in PRIMES}
    good = {p: set() for p in PRIMES}
    degenerate = {p: set() for p in PRIMES}
    first = {}
    pattern_counts = {}

    for i, tup in enumerate(tuples):
        first_p = None
        for p in PRIMES:
            if good_reduction(tup, p):
                good[p].add(i)
            else:
                degenerate[p].add(i)
            ok = a2244_possible_original(tup, p)
            if ok:
                possible[p].add(i)
            else:
                obstructs[p].add(i)
                if first_p is None:
                    first_p = p
            if good_reduction(tup, p) and not ok:
                key = tuple(tuple(row[2]) for row in partition_patterns(tup, p))
                pattern_counts[(p, key)] = pattern_counts.get((p, key), 0) + 1
        first[first_p] = first.get(first_p, 0) + 1

    cover_solutions = minimal_cover(obstructs, universe)
    ff_rows = finite_field_surface_counts([3,5,7,11,13,17,19,23,29,31])

    with output_path.open("w") as out:
        out.write("A(2,2,4,4) local obstruction analysis on K3 surface tuples\n")
        out.write(f"input={input_path}\n")
        out.write(f"tuples={len(tuples)}\n")
        out.write(f"primes={PRIMES}\n\n")

        out.write("Criterion at good reduction:\n")
        out.write("For squared branch residues A,B,C,D, an ordering satisfies A2244 iff the four cross-differences for some 2+2 partition have the same quadratic character.\n")
        out.write("The three partitions tested are 12|34, 13|24, 14|23. Zero character indicates bad/repeated reduction for this criterion.\n\n")

        out.write("Per-prime incidence on rational K3 tuples:\n")
        out.write("p good degenerate possible obstructs\n")
        for p in PRIMES:
            out.write(f"{p} {len(good[p])} {len(degenerate[p])} {len(possible[p])} {len(obstructs[p])}\n")
        out.write("\n")

        out.write("First obstruction distribution in prime order:\n")
        for key in sorted(first, key=lambda x: 10**9 if x is None else x):
            out.write(f"{key} {first[key]}\n")
        out.write("\n")

        out.write("Minimal prime subsets covering all tuple obstructions:\n")
        if cover_solutions:
            for combo in cover_solutions[:20]:
                out.write(" ".join(map(str, combo)) + "\n")
            if len(cover_solutions) > 20:
                out.write(f"... {len(cover_solutions)} total minimal subsets of size {len(cover_solutions[0])}\n")
        else:
            out.write("none\n")
        out.write("\n")

        out.write("Finite-field S(F_p) nondegenerate A2244-lift density:\n")
        out.write("p signed_total signed_good signed_good_possible classes_total classes_good classes_good_possible\n")
        for row in ff_rows:
            out.write("{p} {signed_total} {signed_good} {signed_good_possible} {classes_total} {classes_good} {classes_good_possible}\n".format(**row))
        out.write("\n")

        out.write("Representative good-reduction obstruction patterns, grouped by prime.\n")
        out.write("Format: p count patterns_for_partitions_12|34_13|24_14|23. Characters are -1 nonsquare, 1 square.\n")
        by_prime = {}
        for (p, key), count in pattern_counts.items():
            by_prime.setdefault(p, []).append((count, key))
        for p in sorted(by_prime):
            out.write(f"p={p}\n")
            for count, key in sorted(by_prime[p], reverse=True)[:8]:
                out.write(f"  {count} {key}\n")

    print(f"wrote {output_path}")


if __name__ == "__main__":
    main()
