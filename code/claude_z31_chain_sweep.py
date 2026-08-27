#!/usr/bin/env python3
# claude_z31_chain_sweep.py -- rigorous Eisenstein/Alessandri-Coppola torsion-order
# sweep for dim-2 weight-2 trivial-character newforms from stored Hecke TRACES only.
#
# For a dim-2 form with (totally real quadratic) Hecke field and good primes p,q:
#   t_n   = Tr(a_n)                          (LMFDB mf_newforms.traces, 1-indexed)
#   n_p   = Nm(a_p) = (t_p^2 - t_{p^2} - 2p)/2            [a_{p^2} = a_p^2 - p]
#   z_p   = 2 a_p - t_p  (trace-0 part),  d_p := z_p^2 = t_p^2 - 4 n_p  (integer >= 0)
#   c_pq  = z_p z_q = 2 t_{pq} - t_p t_q   (p != q good)  [a_p a_q = a_{pq}]
#   => for q with sq <= len(traces), seed s with d_s > 0:
#      z_q^2 = c_sq^2 / d_s  (must be a nonnegative integer),
#      n_q   = (t_q^2 - z_q^2)/4  (integrality checked; cross-seed consistency checked)
#   #A_f(F_p) = Nm(1 + p - a_p) = (1+p)^2 - (1+p) t_p + n_p
#   M_odd := odd part of gcd over good odd p of #A_f(F_p) -- a rigorous multiple of
#   the odd torsion order of every member of the isogeny class of A_f.
#
# Validation mode: run on the level<=10^4 exact-tier forms (traces truncated to 100
# to mimic the high-level data) and require exactM | chainM for every form, where
# exactM is the ap-vector result of the 2026-07-30 sweep (mf_sweep_all.csv).
import sys, re, csv, math
from math import gcd
csv.field_size_limit(sys.maxsize)

def parse_nums(s):
    return [int(x) for x in re.findall(r"-?\d+", s)]

def sqfree_core(n):
    n = abs(n); c = 1; d = 2
    while d*d <= n:
        e = 0
        while n % d == 0: n //= d; e += 1
        if e % 2: c *= d
        d += 1
    return c*n if n > 1 else c

def chain_norms(N, t, L):
    """Return dict p -> n_p for good primes p, from traces t (1-indexed list)."""
    P = [p for p in (2,3,5,7,11,13,17,19,23,29,31,37,41,43,47) if N % p != 0]
    tt = lambda n: t[n-1]
    npexact = {}
    dp = {}
    for p in P:
        if p*p <= L:
            num = tt(p)**2 - tt(p*p) - 2*p
            if num % 2: return None, "parity"
            npexact[p] = num//2
            dp[p] = tt(p)**2 - 4*npexact[p]
    seeds = [s for s in npexact if dp[s] > 0]
    for q in P:
        if q in npexact: continue
        vals = set()
        for s in seeds:
            if s == q or s*q > L: continue
            c = 2*tt(s*q) - tt(s)*tt(q)
            if (c*c) % dp[s]: continue  # inconsistent seed -> skip (flag via vals empty)
            zq2 = (c*c)//dp[s]
            num = tt(q)**2 - zq2
            if num < 0 or num % 4: continue
            vals.add(num//4)
        if len(vals) == 1:
            npexact[q] = vals.pop()
        elif len(vals) > 1:
            return None, f"seed-conflict p={q}"
    return npexact, None

def M_odd_of(N, t):
    L = len(t)
    norms, err = chain_norms(N, t, L)
    if norms is None: return None, 0, err
    M = 0; used = 0
    for p, n_p in sorted(norms.items()):
        if p == 2: continue
        loc = (1+p)**2 - (1+p)*t[p-1] + n_p
        M = gcd(M, loc); used += 1
    while M and M % 2 == 0: M //= 2
    return M, used, None

def run(infile, outfile):
    out = csv.writer(open(outfile, "w"))
    out.writerow(["label","level","field_disc","M_odd","n_primes","flag"])
    n_done = n_flag = 0
    for row in csv.DictReader(open(infile)):
        label = row["label"]; N = int(row["level"])
        fp = parse_nums(row["field_poly"])       # ascending c0,c1,c2
        disc = fp[1]*fp[1] - 4*fp[2]*fp[0]
        t = parse_nums(row["traces"])
        M, used, err = M_odd_of(N, t)
        out.writerow([label, N, sqfree_core(disc), M if M is not None else "",
                      used, err or ""])
        n_done += 1; n_flag += (err is not None)
    print(f"{infile}: {n_done} forms, {n_flag} flagged")

def validate(lofile, exactfile):
    exact = {}
    for row in csv.DictReader(open(exactfile)):
        if row["method"].startswith("ap_"):
            exact[row["label"]] = int(row["M_f"])
    bad = tested = flagged = 0
    for row in csv.DictReader(open(lofile)):
        label = row["label"]
        if label not in exact: continue
        N = int(row["level"]); t = parse_nums(row["traces"])[:100]
        M, used, err = M_odd_of(N, t)
        if err or M is None: flagged += 1; continue
        tested += 1
        eodd = exact[label]
        while eodd % 2 == 0: eodd //= 2
        if M % eodd != 0 if eodd else False:
            bad += 1
            if bad <= 10: print(f"VALFAIL {label}: exact_odd={eodd} chain={M}")
    print(f"VALIDATION: {tested} tested, {bad} divisibility failures, {flagged} flagged-skip")
    return bad == 0

if __name__ == "__main__":
    mode = sys.argv[1]
    if mode == "validate":
        ok = validate(sys.argv[2], sys.argv[3])
        sys.exit(0 if ok else 1)
    else:
        run(sys.argv[2], sys.argv[3])
