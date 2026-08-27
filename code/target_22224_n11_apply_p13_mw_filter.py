#!/usr/bin/env python3
"""Apply N11's exact modulo-13^2 MW denominator-depth condition.

For the coefficient basis used by target_22224_n11_rank2_mwsieve.m the
local image is Z/2 x Z/156.  Writing (m,n,ti), the t-denominator has
13-adic valuation at least two only if

    m is even,
    n + 146*(m/2) + torsion_log[ti] is 4 mod 6,

where torsion_log=(0,117,78,39).  The full mod-156 list is one residue
class modulo 6, so this condensed test is exact at depth 13^2; the class
154 has higher simultaneous N,D valuation and is conservatively retained.
"""

from __future__ import annotations
import argparse
import csv
from pathlib import Path

TORSION_LOG = {1: 0, 2: 117, 3: 78, 4: 39}


def allowed(m: int, n: int, ti: int) -> bool:
    if m % 2:
        return False
    j = n + 146 * (m // 2) + TORSION_LOG[ti]
    return j % 6 == 4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    args = ap.parse_args()
    inp, out = Path(args.input), Path(args.output)
    rows = []
    with inp.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if allowed(int(r["m"]), int(r["n"]), int(r["torsion_coset"])):
                rows.append(r)
    with out.open("w") as f:
        f.write("m\tn\ttorsion_coset\n")
        for r in rows:
            f.write(f'{r["m"]}\t{r["n"]}\t{r["torsion_coset"]}\n')
    print(f"input={inp} survivors={sum(1 for _ in inp.open())-1} "
          f"p13_depth2={len(rows)} output={out}")


if __name__ == "__main__":
    main()
