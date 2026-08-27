#!/usr/bin/env python3
"""Fast bounded scan of N11's filtered rank-two MW coefficient lattice."""
from __future__ import annotations
import argparse, csv
from collections import defaultdict
from pathlib import Path

USED = {29, 37, 41, 101, 107, 137, 191}


def bounds(residue: int, modulus: int, n: int):
    lo = (-n - residue + modulus - 1) // modulus
    hi = (n - residue) // modulus
    return lo, hi


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bound", type=int, default=100000)
    ap.add_argument("--lattice", default="results/target_22224_n11_rank2_periodic_lattice_filtered2.tsv")
    ap.add_argument("--profiles", default="results/target_22224_n11_rank2_profiles_p199_all.tsv")
    ap.add_argument("--output", default="results/target_22224_n11_rank2_modular_N100000_p199.tsv")
    args = ap.parse_args()

    allowed = defaultdict(set)
    orders = {}
    totals = {}
    with open(args.profiles) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            p = int(r["prime"])
            if p in USED:
                continue
            o1, o2 = int(r["ord_g1"]), int(r["ord_g2"])
            orders[p] = (o1, o2)
            allowed[p].add((int(r["torsion_coset"]), int(r["m"]), int(r["n"])))
    # Strongest profiles first.
    primes = sorted(allowed, key=lambda p: len(allowed[p]) / (4 * orders[p][0] * orders[p][1]))
    lattices = list(csv.DictReader(open(args.lattice), delimiter="\t"))
    survivors = []
    enumerated = 0
    for c in lattices:
        ti = int(c["torsion_coset"])
        rm, mm = int(c["m_residue"]), int(c["m_modulus"])
        rn, mn = int(c["n_residue"]), int(c["n_modulus"])
        am, bm = bounds(rm, mm, args.bound)
        an, bn = bounds(rn, mn, args.bound)
        for im in range(am, bm + 1):
            m = rm + mm * im
            for jn in range(an, bn + 1):
                n = rn + mn * jn
                enumerated += 1
                good = True
                for p in primes:
                    o1, o2 = orders[p]
                    if (ti, m % o1, n % o2) not in allowed[p]:
                        good = False
                        break
                if good:
                    survivors.append((m, n, ti))
    out = Path(args.output)
    with out.open("w") as f:
        f.write("m\tn\ttorsion_coset\n")
        for row in survivors:
            f.write("\t".join(map(str, row)) + "\n")
    print(f"bound={args.bound} lattices={len(lattices)} profiles={len(primes)} "
          f"enumerated={enumerated} survivors={len(survivors)} output={out}")


if __name__ == "__main__":
    main()
