#!/usr/bin/env python3
"""Enhanced Pass B: recover Norm(a_q) for primes q <= 47 from the trace vector
Tr(a_n), n <= 100, of a dim-2 weight-2 trivial-character newform.

Method: write a_p = (t_p + u_p sqrt(D))/2, D = fundamental disc of Hecke field.
 - direct (good p in {2,3,5,7}): Norm(a_p) = (t_p^2 - T_{p^2} - 2p)/2, so
   D*u_p^2 = t_p^2 - 4*Norm(a_p).
 - chain (good q, anchor good p with u_p != 0, pq <= 100):
   2*T_{pq} - t_p*t_q = D*u_p*u_q  ==> u_q, hence Norm(a_q) = (t_q^2 - D u_q^2)/4.
Norms n_p = (1+p)^2 - (1+p)*t_p + Norm(a_p) = #A_f(F_p).
M_odd := odd part of gcd(n_p) over ALL determined good p (p=2 legal for odd torsion).
Also validates the method against exact ap-based norms on an overlap sample.
"""
import json, csv
from math import gcd, isqrt

SCRATCH = "/tmp/claude-1000/-home-claude-torsion-jac/e62605ac-ce24-4f3b-bad6-68ec708d4962/scratchpad"
PRIMES47 = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47]

def squarefree_kernel(n):
    n = abs(n); k = 1; d = 2
    while d*d <= n:
        e = 0
        while n % d == 0: n //= d; e += 1
        if e % 2: k *= d
        d += 1
    return k * n

def fund_disc(b, c):
    Dp = b*b - 4*c
    sgn = -1 if Dp < 0 else 1
    k = squarefree_kernel(Dp) * sgn
    return k if k % 4 == 1 else 4*k

def perfect_sqrt(n):
    if n < 0: return None
    r = isqrt(n)
    return r if r*r == n else None

def analyze(N, D, tr):
    """returns dict p -> (t_p, Norm(a_p)) for determined good primes p<=47, plus flags"""
    out = {}
    u = {}      # p -> |u_p| for anchors (sign fixed positive for first nonzero anchor)
    flags = []
    # direct primes
    for p in (2,3,5,7):
        if N % p == 0: continue
        t = tr[p-1]; t2 = tr[p*p-1]
        num = t*t - t2 - 2*p
        if num % 2 != 0: flags.append(f"parity_fail_p{p}"); continue
        Np = num//2
        x = t*t - 4*Np           # = D u_p^2
        if x % D != 0: flags.append(f"disc_fail_p{p}"); continue
        s = perfect_sqrt(x//D)
        if s is None: flags.append(f"sq_fail_p{p}"); continue
        out[p] = (t, Np); u[p] = s
    # anchor = smallest direct good prime with u!=0 (largest reach)
    anchors = [(p, u[p]) for p in (2,3,5,7) if p in u and u[p] != 0]
    if anchors:
        pa, ua = anchors[0]
        for q in PRIMES47:
            if q in out or N % q == 0 or q == pa or pa*q > 100: continue
            tq = tr[q-1]; tpq = tr[pa*q-1]
            y = 2*tpq - tr[pa-1]*tq          # = D * (±ua) * u_q
            if y % (D*ua) != 0: flags.append(f"chain_fail_q{q}"); continue
            uq = y // (D*ua)                 # sign irrelevant: squared next
            val = tq*tq - D*uq*uq
            if val % 4 != 0: flags.append(f"norm_fail_q{q}"); continue
            out[q] = (tq, val//4)
    else:
        flags.append("no_anchor")
    return out, flags

def M_from_norms(p2norm):
    Modd = 0; Mstd = 0; n_used = 0
    for p, (t, Np) in p2norm.items():
        n = (1+p)**2 - (1+p)*t + Np
        assert n > 0, (p, t, Np)
        n_used += 1
        Modd = gcd(Modd, n)
        if p > 2: Mstd = gcd(Mstd, n)
    while Modd and Modd % 2 == 0: Modd //= 2
    return Modd, Mstd, n_used

def factorize(n):
    n = abs(n); f = {}
    if n <= 1: return f
    d = 2
    while d*d <= n:
        while n % d == 0: f[d] = f.get(d,0)+1; n //= d
        d += 1 if d == 2 else 2
    if n > 1: f[n] = f.get(n,0)+1
    return f

# ---------------- VALIDATION against exact ap data ----------------
from fractions import Fraction
def exact_norms(rec):
    """from dim2_ap.jsonl record: dict p -> #A(F_p) for good p <= 47"""
    PR46 = [2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61,67,71,73,79,83,89,97,101,
            103,107,109,113,127,131,137,139,149,151,157,163,167,173,179,181,191,193,197,199]
    N = int(rec["level"]); fp = [int(x) for x in rec["field_poly"]]
    c, b = fp[0], fp[1]
    if rec["hecke_ring_power_basis"]:
        basis = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1))]
    else:
        nums = rec["hecke_ring_numerators"]; dens = rec["hecke_ring_denominators"]
        basis = [(Fraction(int(nums[i][0]))/int(dens[i]), Fraction(int(nums[i][1]))/int(dens[i])) for i in range(2)]
    out = {}
    for i, p in enumerate(PR46):
        if p > 47: break
        if N % p == 0: continue
        A, B = rec["ap46"][i]
        uu = A*basis[0][0] + B*basis[1][0]
        vv = A*basis[0][1] + B*basis[1][1]
        s = Fraction(1+p) - uu
        nm = s*s + s*vv*b + vv*vv*c
        assert nm.denominator == 1
        out[p] = int(nm)
    return out

ap_by_label = {}
with open(f"{SCRATCH}/dim2_ap.jsonl") as fh:
    for line in fh:
        r = json.loads(line)
        if int(r["level"]) <= 2500:
            ap_by_label[r["label"]] = r

n_cmp = n_form_ok = n_form_mismatch = n_norms_cmp = 0
mismatches = []
with open(f"{SCRATCH}/valid_sample.jsonl") as fh:
    for line in fh:
        r = json.loads(line)
        lab = r["label"]
        if lab not in ap_by_label: continue
        N = int(r["level"]); fp = [int(x) for x in r["field_poly"]]
        D = fund_disc(fp[1], fp[0])
        tr = [int(x) for x in r["tr100"]]
        got, flags = analyze(N, D, tr)
        exact = exact_norms(ap_by_label[lab])
        ok = True
        for p, (t, Np) in got.items():
            n = (1+p)**2 - (1+p)*t + Np
            n_norms_cmp += 1
            if exact.get(p) != n:
                ok = False; mismatches.append((lab, p, n, exact.get(p)))
        n_cmp += 1
        if ok: n_form_ok += 1
        else: n_form_mismatch += 1
print(f"VALIDATION: {n_cmp} forms compared, {n_norms_cmp} individual norms compared")
print(f"  forms fully consistent: {n_form_ok}, mismatched: {n_form_mismatch}")
for m in mismatches[:20]: print("  MISMATCH:", m)

# ---------------- RUN on traces-only forms ----------------
rows = []
flagcount = {}
with open(f"{SCRATCH}/dim2_traces.jsonl") as fh:
    for line in fh:
        r = json.loads(line)
        lab = r["label"]; N = int(r["level"])
        fp = [int(x) for x in r["field_poly"]]
        D = fund_disc(fp[1], fp[0])
        tr = [int(x) for x in r["tr100"]]
        got, flags = analyze(N, D, tr)
        for fl in flags: flagcount[fl] = flagcount.get(fl,0)+1
        Modd, Mstd, n_used = M_from_norms(got)
        rows.append((lab, N, D, Modd, Mstd, n_used, ";".join(flags)))

print(f"\nPass B enhanced: {len(rows)} forms")
print("flag counts:", flagcount)
from collections import Counter
cnt = Counter()
for lab, N, D, Modd, Mstd, n_used, fl in rows:
    f = factorize(Modd)
    lpf = max([p for p in f if p > 2], default=1)
    cnt[lpf] += 1
print("odd-lpf histogram:", dict(sorted(cnt.items())))
nu = Counter(r[5] for r in rows)
print("n_primes_used histogram:", dict(sorted(nu.items())))

with open(f"{SCRATCH}/pass_b_enhanced.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["label","level","field_disc","M_odd","M_std_oddp","n_primes_used","flags","largest_odd_prime_factor","factorization_Modd"])
    for lab, N, D, Modd, Mstd, n_used, fl in rows:
        f = factorize(Modd)
        lpf = max([p for p in f if p > 2], default=1)
        fstr = "*".join(f"{p}^{e}" if e>1 else str(p) for p,e in sorted(f.items()))
        w.writerow([lab, N, D, Modd, Mstd, n_used, fl, lpf, fstr])
print(f"wrote {SCRATCH}/pass_b_enhanced.csv")

big = sorted(((r) for r in rows if max([p for p in factorize(r[3]) if p>2], default=1) >= 17),
             key=lambda r: (-max([p for p in factorize(r[3]) if p>2], default=1), r[1]))
print(f"\nforms with odd lpf >= 17: {len(big)}")
for lab, N, D, Modd, Mstd, n_used, fl in big[:60]:
    lpf = max([p for p in factorize(Modd) if p>2], default=1)
    print(f"  lpf={lpf:5d} {lab:16s} N={N:7d} D={D:4d} M_odd={Modd:8d} n_p={n_used} {fl}")
