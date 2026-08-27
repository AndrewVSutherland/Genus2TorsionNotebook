#!/usr/bin/env python3
"""Q2b Eisenstein sweep: M_f = gcd over good odd p of Norm(1+p-a_p) for all
dim-2 weight-2 trivial-character newforms in the LMFDB.

Pass A: forms with mf_hecke_nf ap data (p <= 199, exact, handles non-power basis).
Pass B: traces-only forms (level > 10000): p in {3,5,7} via
        Norm(a_p) = (Tr(a_p)^2 - Tr(a_{p^2}) - 2p)/2  (weight 2, trivial char, good p).
"""
import json, sys, csv
from fractions import Fraction
from math import gcd, isqrt

SCRATCH = "/tmp/claude-1000/-home-claude-torsion-jac/e62605ac-ce24-4f3b-bad6-68ec708d4962/scratchpad"

def primes_upto(n):
    s = list(range(n+1)); s[0]=s[1]=0
    for i in range(2, isqrt(n)+1):
        if s[i]:
            for j in range(i*i, n+1, i): s[j]=0
    return [p for p in s if p]

PRIMES46 = primes_upto(199)
assert len(PRIMES46) == 46

def squarefree_kernel(n):
    n = abs(n); k = 1; d = 2
    while d*d <= n:
        e = 0
        while n % d == 0: n //= d; e += 1
        if e % 2: k *= d
        d += 1
    return k * n  # n now 1 or leftover prime

def fund_disc(b, c):
    """fundamental discriminant of Q[x]/(x^2+bx+c)"""
    D = b*b - 4*c
    sgn = -1 if D < 0 else 1
    k = squarefree_kernel(D) * sgn
    return k if k % 4 == 1 else 4*k

def factorize(n):
    n = abs(n); f = {}
    if n <= 1: return f
    d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d,0)+1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n,0)+1
    return f

def norm_quadratic(s, v, b, c):
    """Norm(s - v*beta) where beta root of x^2+bx+c: s^2 + s*v*b + v^2*c"""
    return s*s + s*v*b + v*v*c

rowsA = []
nbad = 0
with open(f"{SCRATCH}/dim2_ap.jsonl") as fh:
    for line in fh:
        r = json.loads(line)
        label = r["label"]; N = int(r["level"])
        fp = [int(x) for x in r["field_poly"]]
        assert len(fp) == 3 and fp[2] == 1, (label, fp)
        c, b = fp[0], fp[1]          # x^2 + b x + c
        ap = r["ap46"]
        assert len(ap) == 46, (label, len(ap))
        if r["hecke_ring_power_basis"]:
            basis = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))]
        else:
            nums = r["hecke_ring_numerators"]; dens = r["hecke_ring_denominators"]
            basis = []
            for i in range(2):
                d = Fraction(int(dens[i]))
                basis.append((Fraction(int(nums[i][0]))/d, Fraction(int(nums[i][1]))/d))
        M = 0; used = 0
        for i, p in enumerate(PRIMES46):
            if p == 2 or N % p == 0: continue
            A, B = ap[i]
            u = A*basis[0][0] + B*basis[1][0]
            v = A*basis[0][1] + B*basis[1][1]
            s = Fraction(1 + p) - u
            nm = norm_quadratic(s, v, b, c)
            assert nm.denominator == 1, (label, p, nm)
            nm = int(nm)
            assert nm > 0, (label, p, nm)
            M = gcd(M, nm); used += 1
            if M == 1: break
        rowsA.append((label, N, fund_disc(b, c), M, used, "ap_p199"))

rowsB = []
with open(f"{SCRATCH}/dim2_traces.jsonl") as fh:
    for line in fh:
        r = json.loads(line)
        label = r["label"]; N = int(r["level"])
        fp = [int(x) for x in r["field_poly"]]
        assert len(fp) == 3 and fp[2] == 1, (label, fp)
        c, b = fp[0], fp[1]
        tr = [int(x) for x in r["tr100"]]
        assert len(tr) == 100 and tr[0] == 2, (label, len(tr), tr[0])
        M = 0; used = 0
        for p in (3, 5, 7):
            if N % p == 0: continue
            Tp = tr[p-1]; Tp2 = tr[p*p-1]
            num = Tp*Tp - Tp2 - 2*p
            assert num % 2 == 0, (label, p)
            Np = num // 2
            nm = (1+p)**2 - (1+p)*Tp + Np
            assert nm > 0, (label, p, nm)
            M = gcd(M, nm); used += 1
        rowsB.append((label, N, fund_disc(b, c), M, used, f"traces_p7_n{used}"))

allrows = rowsA + rowsB
out = f"{SCRATCH}/mf_sweep_all.csv"
with open(out, "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["label","level","field_disc","M_f","n_primes_used","method","largest_prime_factor","factorization"])
    for label, N, D, M, used, meth in allrows:
        f = factorize(M)
        lpf = max(f) if f else (0 if M == 0 else 1)
        fstr = "*".join(f"{p}^{e}" if e > 1 else str(p) for p, e in sorted(f.items())) if f else str(M)
        w.writerow([label, N, D, M, used, meth, lpf, fstr])

# quick console summary
def report(rows, name):
    from collections import Counter
    cnt = Counter()
    big = []
    for label, N, D, M, used, meth in rows:
        f = factorize(M)
        lpf = max(f) if f else 1
        cnt[lpf] += 1
        if lpf >= 17: big.append((lpf, label, N, D, M, used))
    print(f"== {name}: {len(rows)} forms")
    print("largest-prime-factor histogram (lpf: count):")
    for k in sorted(cnt): print(f"  {k}: {cnt[k]}")
    big.sort(reverse=True)
    print(f"forms with lpf >= 17: {len(big)}")
    for lpf, label, N, D, M, used in big[:80]:
        print(f"  lpf={lpf:5d} {label:16s} N={N:7d} disc={D:6d} M_f={M} (n_p={used})")
    return big

bigA = report(rowsA, "Pass A (ap, p<=199)")
bigB = report(rowsB, "Pass B (traces, p in 3,5,7)")

# sanity checks
d = {r[0]: r for r in allrows}
print("\nSanity: 1830.2.a.q ->", d.get("1830.2.a.q"))
print("Sanity: 23.2.a.a   ->", d.get("23.2.a.a"))
print("Sanity: 67.2.a.b   ->", d.get("67.2.a.b"))
print(f"\nwrote {out}")
