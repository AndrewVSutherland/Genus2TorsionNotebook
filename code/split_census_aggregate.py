#!/usr/bin/env python3
"""Aggregate split-census sources into a draft longtable."""
import re, glob, json, subprocess

# 1. bielliptic sweep hits: best (smallest max|coeff|) witness per group
hits = {}
for fn in glob.glob('results/split_bielliptic_H*_s*.log'):
    for line in open(fn):
        m = re.match(r'HIT \[(-?\d+),(-?\d+),(-?\d+),(-?\d+)\] TORSION \[([ 0-9,]*)\] order (\d+)', line)
        if not m: continue
        a,b,c,d = map(int, m.groups()[:4])
        grp = '[' + m.group(5).replace(' ','') + ']'
        ht = max(abs(a),abs(b),abs(c),abs(d))
        if grp not in hits or ht < hits[grp][0]:
            hits[grp] = (ht, (a,b,c,d))
print("bielliptic sweep groups:")
for g in sorted(hits, key=lambda s: (len(eval(s)), eval(s))):
    ht,(a,b,c,d) = hits[g]
    print(f"  {g}  y^2 = {a}x^6 + {b}x^4 + {c}x^2 + {d}")
