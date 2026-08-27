#!/usr/bin/env python3
"""Apply the V1 p=13 depth-five coefficient mask to an MW TSV.

Input rows use the matching basis G1=(3603,216600), G2=(93,-2100).
For m=r+4a and n=s+4b, the local scalar is

    q = a + 10123*b  (mod 13^4).
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

PERIOD = 13**4
CFORMAL = 10123


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input",
                    default="results/target_22224_V1_rank2_modular_N10000.tsv")
    ap.add_argument("--mask",
                    default="results/target_22224_V1_p13_coefficient_mask.tsv")
    ap.add_argument("--output",
                    default="results/target_22224_V1_rank2_modular_N10000_p13deep.tsv")
    args = ap.parse_args()

    with Path(args.mask).open() as f:
        mask = {
            (int(r["m_mod4"]), int(r["n_mod4"]), int(r["q_residue"]))
            for r in csv.DictReader(f, delimiter="\t")
        }

    kept = []
    total = 0
    with Path(args.input).open() as f:
        for row in csv.DictReader(f, delimiter="\t"):
            total += 1
            m, n = int(row["m"]), int(row["n"])
            r, s = m % 4, n % 4
            a, b = (m - r) // 4, (n - s) // 4
            q = (a + CFORMAL * b) % PERIOD
            if (r, s, q) in mask:
                kept.append((m, n, int(row["torsion_coset"]), r, s, q))

    out = Path(args.output)
    with out.open("w") as f:
        f.write("m\tn\ttorsion_coset\tm_mod4\tn_mod4\tq_mod13pow4\n")
        for row in kept:
            f.write("\t".join(map(str, row)) + "\n")
    print(
        f"V1_P13_COEFFICIENT_FILTER input={args.input} total={total} "
        f"kept={len(kept)} output={out}"
    )


if __name__ == "__main__":
    main()
