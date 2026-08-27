#!/usr/bin/env python3
"""Write live ordered p=13 residue bases for M(2,2,2,4)+3.

A residue tuple (a,b,c,d) is kept if it has a cover-open mod-13 cubic-contact
solution which lifts to 13^2 and can make every boundary label exact at first
order.  This is the residue-level input for m2224_plus3_residue_enum.m.
"""

from __future__ import annotations

import argparse
import sys
from itertools import product
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from m2224_plus3_boundary13_lift import (  # noqa: E402
    P,
    boundary_labels,
    can_make_all_labels_first_order,
    coeffs_mod,
    cover_open,
    label_key,
    linear_lift_space,
    raw_contact_witnesses_by_key,
)


def live_residues() -> dict[tuple[int, int, int, int], set[str]]:
    by_key = raw_contact_witnesses_by_key()
    live: dict[tuple[int, int, int, int], set[str]] = {}
    for vals in product(range(P), repeat=4):
        labels = boundary_labels(vals)
        if labels == ("open",):
            continue
        key = coeffs_mod(*vals)
        for wit in by_key.get(key, []):
            if not cover_open(key, wit):
                continue
            z = (*vals, *wit)
            space = linear_lift_space(z)
            if space is None:
                continue
            y0, basis = space
            if can_make_all_labels_first_order(labels, y0, basis):
                live.setdefault(tuple(vals), set()).add(label_key(labels))
    return live


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", default="data/m2224_plus3_live13_residues.txt")
    args = parser.parse_args()

    live = live_residues()
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as f:
        f.write("# a,b,c,d | signatures\n")
        for vals in sorted(live):
            sigs = ";".join(sorted(live[vals]))
            f.write(f"{vals[0]},{vals[1]},{vals[2]},{vals[3]} | {sigs}\n")
    print("live_residue_bases", len(live), "out", out)


if __name__ == "__main__":
    main()
