#!/usr/bin/env python3
"""
claude_ov_m612ap_satcompare.py    lane 9 ([6,12])   2026-07-25

Put the two sides of the Prym identification next to each other:

  * results/ap3/ap3_q*.txt          — from the genus-4 curve E8:
        PRYMORD = #J(E8)(F_q) / #J(E4)(F_q) = #Pr(F_q)
        WORDERS = ord(W mod q), W the intrinsic Abel-Prym class
  * results/claude_ov_m612ap_satcheck.log — from the genus-2 curve
        D : y^2 = -3x^6+24x^3-75 :  #J(D)(F_q)  and  ord(g0 mod q),
        g0 = (x^2+2x+4, 5x+5) the generator of J(D)(Q) (claude_ov_m612ap_satindex.m).

Agreement of #Pr(F_q) with #J(D)(F_q) is the point-count half of the isogeny
certificate; agreement of ord(W mod q) with ord(g0 mod q) is evidence that W
corresponds to +-g0, i.e. that W is a GENERATOR of Pr(Q) modulo torsion --
the sieve's saturation hypothesis.

Usage: python3 code/claude_ov_m612ap_satcompare.py
"""
import glob, sys

satlog = "results/claude_ov_m612ap_satcheck.log"
apdir = "results/ap3"

sat = {}
for line in open(satlog):
    t = line.split()
    if len(t) == 4 and t[0].isdigit():
        sat[int(t[0])] = (int(t[1]), int(t[2]))

print("%-5s %-14s %-14s %-14s %-12s %s"
      % ("q", "#J(D)(F_q)", "#Pr(F_q)", "ord(g0 mod q)", "ord(W mod q)", "agree"))
ok, bad = 0, []
for f in sorted(glob.glob(apdir + "/ap3_q*.txt"),
                key=lambda f: int(f.split("_q")[-1].split(".")[0])):
    q = int(f.split("_q")[-1].split(".")[0])
    D = dict(l.strip().split(" ", 1) for l in open(f) if " " in l.strip())
    if "WORDERS" not in D or q not in sat:
        continue
    m = max(int(u) for u in D["WORDERS"].rstrip(",").split(","))
    pr = D.get("PRYMORD", "?")
    hD, og = sat[q]
    good = (og == m) and (str(hD) == pr)
    ok += 1 if good else 0
    if not good:
        bad.append(q)
    print("%-5d %-14d %-14s %-14d %-12d %s"
          % (q, hD, pr, og, m, "YES" if good else "NO"))
print("\nboth quantities agree at %d primes; disagreements: %s" % (ok, bad or "none"))
