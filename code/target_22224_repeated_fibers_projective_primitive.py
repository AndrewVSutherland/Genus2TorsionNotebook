#!/usr/bin/env python3
"""Primitive projective p=11,13 contact masks for P8 and F1--F7.

The family parameter is T=r/s and the signed branch tuple is

    (f1*s^2, f2*s^2, f3*s^2, d0*r^2).

Before reduction we divide by its common local p-power.  This resolves the
two all-zero raw charts (F2 at p=11 infinity and F3 at p=13 infinity) and
P8's hidden p=11 finite T=0 chart.
"""
from __future__ import annotations
import csv
from pathlib import Path

ROOT = Path("results")
FIBERS = {
    "P8": ((528, -726, -891), 50),
    "F1": ((-1470, -630, 336), 25),
    "F2": ((-720, 20, 300), -363),
    "F3": ((-612, 34, 289), -338),
    "F4": ((-126, 28, 49), -50),
    "F5": ((-112, 14, 49), -50),
    "F6": ((-50, 30, 45), -48),
    "F7": ((-18, 1, 16), -50),
}


def vp(x: int, p: int) -> int:
    if x == 0:
        return 100
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


def primitive(vals, p):
    e = min(vp(x, p) for x in vals)
    assert e < 100
    return tuple((x // p**e) % p for x in vals), e


def mask(p):
    ans = set()
    path = ROOT / f"target_22224_direct_contact_deep13_boundary_p{p}.tsv"
    with path.open() as f:
        for r in csv.DictReader(f, delimiter="\t"):
            a = [int(r[k]) % p for k in "abcd"]
            for u in range(1, p):
                ans.add(tuple(sorted((u * x % p) ** 2 % p for x in a)))
    return ans


def contact(vals, p, m):
    key = tuple(sorted((x % p) ** 2 % p for x in vals))
    return key in m


def main():
    rows = []
    lines = []
    for p in (11, 13):
        m = mask(p)
        lines.append(f"p={p} projected_contact_square_multisets={len(m)}")
        for name, (fixed, d0) in FIBERS.items():
            finite = []
            finite_primitive = {}
            for t in range(p):
                vals, e = primitive(tuple(fixed) + (d0 * t * t,), p)
                if contact(vals, p, m):
                    finite.append(t)
                finite_primitive[t] = (vals, e)
            # Infinity disk: r=1, s=p*z with z a unit.  At depth >=2 the
            # first stable representative is s=p^2.
            depth1 = []
            depth1_data = {}
            for z in range(1, p):
                vals, e = primitive(tuple(x * (p * z) ** 2 for x in fixed) + (d0,), p)
                if contact(vals, p, m):
                    depth1.append(z)
                depth1_data[z] = (vals, e)
            vals2, e2 = primitive(tuple(x * p**4 for x in fixed) + (d0,), p)
            depth2 = contact(vals2, p, m)
            rows.append((name, p, ",".join(map(str, finite)),
                         ",".join(map(str, depth1)), int(depth2)))
            exceptional = []
            if finite_primitive[0][1]:
                exceptional.append(f"finite0={finite_primitive[0]}")
            if len(set(v[0] for v in depth1_data.values())) > 1 or e2:
                exceptional.append(f"inf1_z1={depth1_data[1]} inf2={vals2,e2}")
            lines.append(
                f"{name}: finite={finite}; infinity_depth1_z={depth1}; "
                f"infinity_depth>=2={depth2}" +
                (("; " + "; ".join(exceptional)) if exceptional else "")
            )
    out = ROOT / "target_22224_repeated_fibers_projective_primitive.tsv"
    with out.open("w") as f:
        f.write("fiber\tprime\tfinite_T_classes\tinfinity_depth1_z\tinfinity_depth_ge2\n")
        for row in rows:
            f.write("\t".join(map(str, row)) + "\n")
    log = ROOT / "target_22224_repeated_fibers_projective_primitive.log"
    log.write_text("\n".join(lines) + "\n")
    print(log.read_text(), end="")


if __name__ == "__main__":
    main()
