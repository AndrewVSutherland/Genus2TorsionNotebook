#!/usr/bin/env python3
"""
claude_ov_m612ap_regress.py   lane 9 ([6,12])   2026-07-25

Regression check between the two independent sieve workers:
  results/ap2/ap2_q*.txt   from code/claude_ov_m612ap_sieve2.m (previous session)
  results/ap3/ap3_q*.txt   from code/claude_ov_m612ap_sieve3.m (rewrite)
For every prime computed by both, JORDER, WORDERS and RES2 must be identical.
sieve3 changed the generator-skip logic, added RES3/RES9 and fixed a parse
error, so this checks that nothing else moved.

Usage: python3 code/claude_ov_m612ap_regress.py
"""
import glob

def parse(f):
    return dict(l.strip().split(" ", 1) for l in open(f) if " " in l.strip())

same, bad = 0, []
for f in sorted(glob.glob("results/ap2/ap2_q*.txt"),
                key=lambda f: int(f.split("_q")[-1].split(".")[0])):
    q = int(f.split("_q")[-1].split(".")[0])
    try:
        A, B = parse(f), parse("results/ap3/ap3_q%d.txt" % q)
    except FileNotFoundError:
        continue
    if "RES2" not in A or "RES2" not in B:
        continue
    keys = ("JORDER", "JINV", "WORDERS", "SORDER", "NONFIXED", "INS", "RES2",
            "MINMULT", "JE4ORD", "PRYMORD", "DIVISIBLE", "TESTED")
    ok = all(A.get(k) == B.get(k) for k in keys)
    same += ok
    if not ok:
        bad.append((q, [(k, A.get(k), B.get(k)) for k in keys if A.get(k) != B.get(k)]))
print("primes computed by BOTH workers with identical output: %d" % same)
print("differing: %d" % len(bad))
for q, d in bad:
    print("  q = %d : %s" % (q, d))
