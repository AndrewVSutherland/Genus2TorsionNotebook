#!/usr/bin/env python3
"""Exact low-memory search on the split-quadratic-B slice d=5/3.

The Hermite quotient S has a rational root u coming from the linear M(12)
factor. Parametrize (u,v) on v^2=Q(u) by the slope k through (0,s):

    u=(t^2-s^2-1-2*s*k)/(k^2-1),  v=s+k*u.

Then a is recovered from A(u)=v^5. If rho is the root of the M(12) linear
factor, comparison of first Taylor coefficients gives

    S'(u)/lc(S) = (rho*(rho+1)/d^2)^2.

For fixed d=5/3 this is two successive rational-square tests. Remaining
Taylor coefficients recover w, and every hit is checked against the complete
rational polynomial identity S=tau*F(d*Z-c), with tau a nonzero square.

The local hash, intrinsic-coordinate CRT (including leave-one-prime-out), and
optional rational-height sieve all have small, predictable memory use.
"""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import math
import sys
from fractions import Fraction
from pathlib import Path


D_FIXED = Fraction(5, 3)
HERE = Path(__file__).resolve().parent


def import_file(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, HERE / filename)
    assert spec and spec.loader
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


FQ = import_file("fullquad_d53_base", "m12_general5_fullquad_rootpair.py")
C2 = import_file("fullquad_d53_exact", "m12_general5_fullquad_c2_crt.py")


def peval(f, x):
    ans = 0
    for y in reversed(f):
        ans = ans*x+y
    return ans


def taylor_coefficient(f, u, j: int):
    return sum(math.comb(i, j)*f[i]*u**(i-j) for i in range(j, len(f)))


def sqrt_fraction(x: Fraction):
    if x < 0:
        return None
    rn, rd = math.isqrt(x.numerator), math.isqrt(x.denominator)
    if rn*rn != x.numerator or rd*rd != x.denominator:
        return None
    return Fraction(rn, rd)


def compact_exact(b, w, c, d):
    x = [-c, d]
    L = C2.padd([b], C2.pscale(x, 2*b-1))
    H = C2.padd(x, C2.pscale(C2.padd([1], C2.pscale(x, b)), w))
    inner = C2.padd(
        C2.pmul(L, C2.pmul(H, H)),
        C2.pscale(
            C2.pmul(
                C2.pmul(C2.padd([1], x), C2.padd([1], x)),
                C2.padd(C2.pscale(L, w), C2.pscale(C2.pmul(x, x), -1)),
            ),
            4*b,
        ),
    )
    return C2.pmul(L, inner)


def exact_match(vals):
    b, w, c, d, s, t, a = vals
    if not b or b == 1 or 2*b == 1 or not w or not d or not s or not t:
        return None
    S, A, Q = C2.hermite_exact(s, t, a)
    G = compact_exact(b, w, c, d)
    S += [Fraction(0)]*(6-len(S))
    G += [Fraction(0)]*(6-len(G))
    if not G[5]:
        return None
    tau = S[5]/G[5]
    if not tau or any(S[i] != tau*G[i] for i in range(6)):
        return None
    if not C2.is_fraction_square(tau):
        return None
    if Q[1]*Q[1]-4*Q[0] == 0:
        return None
    return tau, A, Q, S, G


def local_points(p: int):
    """All smooth/open signed-endpoint points on d=5/3 over F_p."""
    if p in (2, 3, 5):
        return []
    d = 5*pow(3, -1, p) % p
    table, _ = FQ.build_S_hash(p)
    out = []
    for b in range(p):
        for w in range(p):
            if not FQ.base_open(b, w, p):
                continue
            if not FQ.squarefree(FQ.compact_F_transformed(b, w, 0, 1, p), p):
                continue
            for c in range(p):
                G = FQ.compact_F_transformed(b, w, c, d, p)
                key, gpivot = FQ.projective_key(G, p)
                for rec in table.get(key, []):
                    ans = FQ.check_match(b, w, c, d, rec, G, gpivot, p)
                    if ans is None or ans[3] != 5:
                        continue
                    vals = ans[0]
                    out.append((vals[0], vals[1], vals[2], vals[4], vals[5], vals[6]))
    return sorted(set(out))


def intrinsic_mod(point, p: int):
    b, _w, c, s, t, a = point
    d = 5*pow(3, -1, p) % p
    rho = -b*pow((2*b-1) % p, -1, p) % p
    u = (rho+c)*pow(d, -1, p) % p
    if not u:
        return None
    S, A, Q = FQ.hermite_S(s, t, a, p)
    qv = peval(Q, u) % p
    av = peval(A, u) % p
    if not qv or av*av % p != pow(qv, 5, p):
        return None
    v = av*pow(qv, -2, p) % p
    if v*v % p != qv or pow(v, 5, p) != av:
        return None
    k = (v-s)*pow(u, -1, p) % p
    if (k*k-1) % p == 0:
        return None
    expected = (t*t-s*s-1-2*s*k)*pow((k*k-1) % p, -1, p) % p
    if u != expected or peval(S, u) % p:
        raise ArithmeticError("intrinsic local-coordinate control failed")
    return s, t, k


def recover_candidates(s: Fraction, t: Fraction, k: Fraction):
    if not s or not t or k*k == 1:
        return []
    u = (t*t-s*s-1-2*s*k)/(k*k-1)
    if not u or u == 1:
        return []
    v = s+k*u
    _S0, A0, _Q0 = C2.hermite_exact(s, t, Fraction(0))
    a = (v**5-peval(A0, u))/(u*u*(u-1)*(u-1))
    S, A, Q = C2.hermite_exact(s, t, a)
    S += [Fraction(0)]*(6-len(S))
    if peval(S, u) or peval(Q, u) != v*v or peval(A, u) != v**5:
        raise ArithmeticError("exact root parametrization control failed")
    lc = S[5]
    if not lc:
        return []
    R = taylor_coefficient(S, u, 1)/lc
    h0 = sqrt_fraction(R)
    if h0 is None:
        return []

    hits = []
    for h in ([h0] if not h0 else [h0, -h0]):
        m0 = sqrt_fraction(1+4*h*D_FIXED*D_FIXED)
        if m0 is None:
            continue
        for m in ([m0] if not m0 else [m0, -m0]):
            rho = (m-1)/2
            if not rho or rho == -1 or 2*rho+1 == 0:
                continue
            if rho*(rho+1) != h*D_FIXED*D_FIXED:
                raise ArithmeticError("fixed-d derivative cover failed")
            b = rho/(2*rho+1)
            c = D_FIXED*u-rho

            C2t = taylor_coefficient(S, u, 2)
            C3t = taylor_coefficient(S, u, 3)
            X2 = C2t*D_FIXED**3/lc
            X3 = C3t*D_FIXED**2/lc
            Y2 = 4*rho*(2*rho+1)**2*X2
            Y3 = 2*rho*(2*rho+1)**2*X3
            D23 = (rho+1)**2*Y3-rho*Y2
            K23 = -rho*(2*rho+1)**2*(
                4*rho**4-12*rho**3-30*rho**2-18*rho-3
            )
            qden = (2*rho+1)*(6*rho+1)
            if qden:
                qroot = (D23-K23)/qden
            else:
                C4t = taylor_coefficient(S, u, 4)
                X4 = C4t*D_FIXED/lc
                Y4 = 4*rho*(2*rho+1)**2*X4
                D24 = (rho+1)**4*Y4-rho*rho*Y2
                K24 = (2*rho+1)**5*(6*rho*rho+6*rho+1)
                qden = 6*rho*(2*rho+1)**2
                if not qden:
                    continue
                qroot = (D24-K24)/qden
            w = qroot/(rho+1)**2
            ans = exact_match((b, w, c, D_FIXED, s, t, a))
            if ans is not None:
                row = (b, w, c, D_FIXED, s, t, a, ans[0], u, v, k)
                if row not in hits:
                    hits.append(row)
    return hits


def rational_reconstruction(r: int, modulus: int):
    r %= modulus
    bound = math.isqrt((modulus-1)//2)
    old_r, new_r, old_t, new_t = modulus, r, 0, 1
    while abs(new_r) > bound:
        if not new_r:
            return None
        q = old_r//new_r
        old_r, new_r = new_r, old_r-q*new_r
        old_t, new_t = new_t, old_t-q*new_t
    n, den = new_r, new_t
    if den < 0:
        n, den = -n, -den
    if not den or den > bound or abs(n) > bound or math.gcd(abs(n), den) != 1:
        return None
    if (r*den-n) % modulus:
        return None
    return Fraction(n, den)


def local_intrinsic_rows(primes):
    rows = {}
    for p in primes:
        ambient = local_points(p)
        intrinsic = sorted({q for P in ambient if (q := intrinsic_mod(P, p)) is not None})
        rows[p] = intrinsic
        print(f"LOCAL_D53 p={p} ambient={len(ambient)} unique_stk={len(intrinsic)} data={intrinsic}")
    return rows


def search_crt(primes, rows, allow_one_bad=True):
    if any(not rows[p] for p in primes):
        print("CRT_ABORT locally empty prime")
        return []
    subsets = [tuple(primes)]
    if allow_one_bad:
        subsets += [tuple(q for q in primes if q != p) for p in primes]
    exact_triples, hits = set(), []
    tested = reconstructed = 0
    for subset in subsets:
        modulus = math.prod(subset)
        combos = math.prod(len(rows[p]) for p in subset)
        here_rec = here_new = 0
        cache = {}
        for choice in itertools.product(*(rows[p] for p in subset)):
            tested += 1
            vals = []
            for j in range(3):
                residues = tuple(x[j] for x in choice)
                key = (j, residues)
                if key not in cache:
                    r, M = C2.crt(residues, subset)
                    cache[key] = rational_reconstruction(r, M)
                vals.append(cache[key])
            if any(x is None for x in vals):
                continue
            reconstructed += 1
            here_rec += 1
            triple = tuple(vals)
            if triple in exact_triples:
                continue
            exact_triples.add(triple)
            here_new += 1
            for hit in recover_candidates(*triple):
                if hit not in hits:
                    hits.append(hit)
                    print("EXACT_HIT_D53 b,w,c,d,s,t,a,tau,u,v,k=", hit)
        print(f"CRT_SUBSET primes={subset} modulus={modulus} combinations={combos} "
              f"reconstructed={here_rec} new_exact_triples={here_new}")
    print(f"CRT_D53_DONE subsets={len(subsets)} combinations={tested} "
          f"reconstructed={reconstructed} unique_exact_triples={len(exact_triples)} "
          f"exact_hits={len(hits)}")
    return hits


def rational_parameters(height: int):
    return [Fraction(n, d) for d in range(1, height+1)
            for n in range(-height, height+1) if math.gcd(abs(n), d) == 1]


def search_height(height: int, primes, rows):
    """Exhaust reduced s,t,k of height <=H after a conservative local sieve."""
    params = rational_parameters(height)
    nonzero = {i for i, x in enumerate(params) if x}
    allidx = set(range(len(params)))
    res, bad, byres = {}, {}, {}
    st_by_s, k_by_st = {}, {}
    for p in primes:
        rp, bp, ip = [], set(), {}
        for i, x in enumerate(params):
            if x.denominator % p == 0:
                rp.append(None)
                bp.add(i)
            else:
                y = x.numerator*pow(x.denominator, -1, p) % p
                rp.append(y)
                ip.setdefault(y, set()).add(i)
        res[p], bad[p], byres[p] = rp, bp, ip
        smap, kmap = {}, {}
        for s, t, k in rows[p]:
            smap.setdefault(s, set()).add(t)
            kmap.setdefault((s, t), set()).add(k)
        st_by_s[p], k_by_st[p] = smap, kmap

    checked_pairs = checked_triples = exact_triples = 0
    hits = []
    for is_, s in enumerate(params):
        if not s:
            continue
        tcand = nonzero.copy()
        for p in primes:
            sr = res[p][is_]
            if sr is None:
                continue
            allowed = st_by_s[p].get(sr, set())
            accept = bad[p].copy()
            for tr in allowed:
                accept.update(byres[p].get(tr, ()))
            tcand.intersection_update(accept)
            if not tcand:
                break
        for it in tcand:
            t = params[it]
            checked_pairs += 1
            kcand = allidx.copy()
            for p in primes:
                sr, tr = res[p][is_], res[p][it]
                if sr is None or tr is None:
                    continue
                allowed = k_by_st[p].get((sr, tr), set())
                accept = bad[p].copy()
                for kr in allowed:
                    accept.update(byres[p].get(kr, ()))
                kcand.intersection_update(accept)
                if not kcand:
                    break
            checked_triples += len(kcand)
            for ik in kcand:
                k = params[ik]
                if k*k == 1:
                    continue
                exact_triples += 1
                for hit in recover_candidates(s, t, k):
                    if hit not in hits:
                        hits.append(hit)
                        print("EXACT_HEIGHT_HIT_D53 b,w,c,d,s,t,a,tau,u,v,k=", hit)
        if is_ % 256 == 0:
            print(f"HEIGHT_PROGRESS s_index={is_}/{len(params)} pairs={checked_pairs} "
                  f"triples={checked_triples}")
    formal = len(params)**3
    print(f"HEIGHT_D53_DONE H={height} parameters={len(params)} formal_triples={formal} "
          f"surviving_pairs={checked_pairs} modular_triples={checked_triples} "
          f"exact_tested={exact_triples} exact_hits={len(hits)}")
    return hits


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--primes", default="7,11,13,17,19,23,29")
    ap.add_argument("--mode", choices=["local", "crt", "height", "all"], default="all")
    ap.add_argument("--height", type=int, default=30)
    ap.add_argument("--no-one-bad", action="store_true")
    args = ap.parse_args()
    primes = [int(x) for x in args.primes.split(",") if x]
    rows = local_intrinsic_rows(primes)
    if args.mode in ("crt", "all"):
        search_crt(primes, rows, not args.no_one_bad)
    if args.mode in ("height", "all"):
        search_height(args.height, primes, rows)


if __name__ == "__main__":
    main()
