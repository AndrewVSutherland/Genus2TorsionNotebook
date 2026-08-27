#!/usr/bin/env python3
"""
claude_ov_m612ap_combine.py    lane 9 ([6,12])   2026-07-25

Combine the per-prime output of code/claude_ov_m612ap_sieve2.m into the
Abel-Prym Mordell-Weil sieve conclusion.

For each reference class W_j (j = 1,2,3) and each good prime q the Magma
worker records
    m       = ord(Wbar_j) in J(E8)(F_q),
    RES_j   = { d mod m : some point R of E8(F_q) has [R - iota R] = d*Wbar_j },
where the two iota-fixed boundary points contribute d = 0.

If P is a rational point of E8 with alpha(P) = d * W_j in Pr(Q) (which is
every rational point once <W_j> is saturated in Pr(Q) modulo torsion),
then for every good q the residue of d modulo ord(Wbar_j) lies in RES_j.
Intersecting over q, prime power by prime power, bounds d.

Usage: python3 code/claude_ov_m612ap_combine.py results/ap2 ap2
"""
import sys, os, glob
from math import gcd

def parse(fn):
    D = {}
    for line in open(fn):
        line = line.strip()
        if not line:
            continue
        t = line.split(" ", 1)
        D[t[0]] = t[1] if len(t) > 1 else ""
    return D

def seq(s):
    if s is None or s.strip() in ("", "-"):
        return []
    return [int(u) for u in s.strip().rstrip(",").split(",") if u != ""]

def prime_powers(m):
    out = []
    d = 2
    mm = m
    while d * d <= mm:
        if mm % d == 0:
            e = 0
            while mm % d == 0:
                mm //= d
                e += 1
            out += [(d, d**k) for k in range(1, e + 1)]
        d += 1
    if mm > 1:
        out.append((mm, mm))
    return out

def main():
    wdir = sys.argv[1] if len(sys.argv) > 1 else "results/ap2"
    tag = sys.argv[2] if len(sys.argv) > 2 else "ap2"
    files = sorted(glob.glob(os.path.join(wdir, tag + "_q*.txt")),
                   key=lambda f: int(f.split("_q")[-1].split(".")[0]))
    rows = []
    skipped = []
    for f in files:
        q = int(f.split("_q")[-1].split(".")[0])
        D = parse(f)
        if "SKIP" in D or "JORDER" not in D:
            skipped.append((q, D.get("SKIP", "no JORDER")))
            continue
        rows.append((q, D))

    print("=== Abel-Prym MW sieve combination (E8, [6,12] lane) ===")
    print("usable primes: %d   skipped: %d" % (len(rows), len(skipped)))
    for q, why in skipped:
        print("  SKIP q=%-4d %s" % (q, why))
    if not rows:
        return

    # --- torsion bound and non-torsion certificate ---------------------
    g = 0
    for q, D in rows:
        g = gcd(g, int(D["JORDER"]))
    print("\ngcd over usable q of #J(E8)(F_q)  (bounds #J(E8)(Q)_tors, hence")
    print("  the torsion of Pr(Q)):  %d" % g)

    for j in (1, 2, 3):
        ords = {}
        for q, D in rows:
            w = seq(D.get("WORDERS"))
            if len(w) >= j:
                ords[q] = w[j - 1]
        vals = sorted(set(ords.values()))
        print("\n--- reference W%d ---" % j)
        print("  ord(W%d mod q) takes %d distinct values, e.g. %s"
              % (j, len(vals), vals[:8]))
        if len(vals) > 1:
            qs = [q for q in ords if ords[q] == vals[0]][:1] + \
                 [q for q in ords if ords[q] == vals[-1]][:1]
            print("  => W%d has INFINITE order in Pr(Q): reduction is injective on"
                  % j)
            print("     torsion, but ord(W%d mod %d) = %d != %d = ord(W%d mod %d)."
                  % (j, qs[0], ords[qs[0]], ords[qs[1]], j, qs[1]))
        else:
            print("  => order constant; W%d may be torsion (order dividing %d)."
                  % (j, vals[0]))

        # --- the sieve ----------------------------------------------------
        allowed = {}          # prime power -> set of allowed residues
        users = {}
        for q, D in rows:
            m = ords.get(q, 0)
            if m <= 1:
                continue
            R = set(seq(D.get("RES%d" % j)))
            for (ell, pk) in prime_powers(m):
                A = set(r % pk for r in R)
                if pk in allowed:
                    allowed[pk] &= A
                    users[pk].append(q)
                else:
                    allowed[pk] = A
                    users[pk] = [q]
        M = 1
        best = {}
        for pk in sorted(allowed):
            A = allowed[pk]
            ell = prime_powers(pk)[0][0]
            print("    modulus %-8d  %2d primes  |allowed| = %-5d %s"
                  % (pk, len(users[pk]), len(A),
                     sorted(A) if len(A) <= 10 else ""))
            if A == {0}:
                if ell not in best or pk > best[ell]:
                    best[ell] = pk
        for ell in sorted(best):
            M *= best[ell]
        print("  SIEVE CONCLUSION (reference W%d): every rational point P of E8" % j)
        print("    with alpha(P) = d*W%d satisfies  d = 0 mod %d" % (j, M))
        if M > 1:
            print("    (product of the prime powers whose allowed set collapsed to {0})")

if __name__ == "__main__":
    main()
