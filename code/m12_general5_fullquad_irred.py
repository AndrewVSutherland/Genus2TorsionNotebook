#!/usr/bin/env python3
"""Low-memory nonsplit/full quadratic-B census for M(12)+5.

This is the quadratic-algebra companion to ``m12_general5_fullquad_rootpair``.
It does not choose roots of B.  Normalize

    B = x^2+r*x+n,        u=s*x+t,
    q = u^2+(1-s^2)*B.

If gcd(u,B)=1, Hensel lifting gives the unique square root ``v`` of q modulo
``B^2`` which reduces to u modulo B.  Every monic quintic A satisfying
``A^2=q^5 (mod B^2)`` is then

    A = (q^2*v mod B^2) + B^2*(x+a).

Consequently ``S=(A^2-q^5)/B^2`` is a quintic.  The desired norm identity is
equivalent to

    S = tau*F_{b,w},       tau a nonzero square.

Unlike the root-pair chart, this parametrization covers irreducible quadratic
B over Q (and nonsplit B over finite fields).  The finite-field routines hash
only the O(p^2) M(12) quintics and stream the O(p^5) quadratic-algebra side,
so memory stays small.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from typing import Iterable, Sequence

from m12_general5_fullquad_rootpair import (
    add,
    base_open,
    compact_F_transformed,
    divmod_poly,
    gcd_poly,
    is_square_nonzero,
    mul,
    pad,
    power,
    projective_key,
    scale,
    squarefree,
    sub,
)


def legendre(a: int, p: int) -> int:
    a %= p
    if not a:
        return 0
    z = pow(a, (p - 1) // 2, p)
    return -1 if z == p - 1 else z


def quadratic_data(r: int, n: int, s: int, t: int, a: int, p: int):
    """Return (S,A,q,B), or None off the coprime/squarefree chart."""
    r %= p; n %= p; s %= p; t %= p; a %= p
    B = [n, r, 1]
    B2 = mul(B, B, p)
    k = (1 - s*s) % p
    u = [t, s]

    # Norm_{F_p[x]/B}(s*x+t).  It is also Res(B,u).
    D = (t*t - r*s*t + n*s*s) % p
    if D == 0:
        return None

    # q=u^2+kB is monic.  h=(k/2)u^{-1} mod B and v=u+B*h
    # satisfy v^2=q mod B^2.
    q = add(mul(u, u, p), scale(B, k, p), p)
    q = list(pad(q, 3))
    if not squarefree(q, p):
        return None
    fac = k * pow(2*D % p, -1, p) % p
    hinv = [fac*(t-r*s) % p, -fac*s % p]
    v = add(u, mul(B, hinv, p), p)
    vv_minus_q = sub(mul(v, v, p), q, p)
    _, rem = divmod_poly(vv_minus_q, B2, p)
    assert rem == [0]

    q2v = mul(mul(q, q, p), v, p)
    A0 = divmod_poly(q2v, B2, p)[1]
    A = add(A0, mul(B2, [a, 1], p), p)
    A = list(pad(A, 6))
    numerator = sub(mul(A, A, p), power(q, 5, p), p)
    S, rem = divmod_poly(numerator, B2, p)
    assert rem == [0]
    S = list(pad(S, 6))
    return S, A, q, B


def build_F_hash(p: int):
    table = defaultdict(list)
    curves = 0
    for b in range(p):
        for w in range(p):
            if not base_open(b, w, p):
                continue
            F = compact_F_transformed(b, w, 0, 1, p)
            if not squarefree(F, p):
                continue
            curves += 1
            key, pivot = projective_key(F, p)
            table[key].append((b, w, pivot, F))
    return table, curves


def match_records(r: int, n: int, s: int, t: int, a: int, p: int,
                  Ftable):
    data = quadratic_data(r, n, s, t, a, p)
    if data is None:
        return []
    S, A, q, B = data
    if S[5] == 0 or not squarefree(S, p):
        return []
    key, _ = projective_key(S, p)
    ans = []
    for b, w, _, F in Ftable.get(key, ()):
        tau = S[5] * pow(F[5], -1, p) % p
        if not is_square_nonzero(tau, p):
            continue
        if len(gcd_poly(q, F, p)) != 1:
            continue
        ans.append((b, w, r % p, n % p, s % p, t % p, a % p,
                    tau, q, B, A, S, F))
    return ans


def full_census(p: int, discriminant: str = "nonsplit",
                sample_limit: int = 4):
    Ftable, curves = build_F_hash(p)
    tested_B = admissible = matches = 0
    masks = set()
    samples = []
    for r in range(p):
        for n in range(p):
            chi = legendre(r*r-4*n, p)
            if discriminant == "nonsplit" and chi != -1:
                continue
            if discriminant == "split" and chi != 1:
                continue
            if discriminant == "separable" and chi == 0:
                continue
            tested_B += 1
            hit_B = False
            for s in range(p):
                for t in range(p):
                    # Build once with a=0 to count admissible (r,n,s,t).
                    if quadratic_data(r, n, s, t, 0, p) is None:
                        continue
                    admissible += 1
                    for a in range(p):
                        recs = match_records(r, n, s, t, a, p, Ftable)
                        if not recs:
                            continue
                        matches += len(recs)
                        hit_B = True
                        if len(samples) < sample_limit:
                            samples.extend(recs[:sample_limit-len(samples)])
            if hit_B:
                masks.add((r, n))
    return dict(p=p, discriminant=discriminant,
                smooth_base_curves=curves, tested_B=tested_B,
                admissible_Bst=admissible, open_matches=matches,
                B_mask=sorted(masks), samples=samples)


def fixed_B_census(p: int, r: int, n: int, sample_limit: int = 4,
                   Ftable=None):
    if Ftable is None:
        Ftable, curves = build_F_hash(p)
    else:
        curves = None
    matches = 0
    points = []
    for s in range(p):
        for t in range(p):
            if quadratic_data(r, n, s, t, 0, p) is None:
                continue
            for a in range(p):
                recs = match_records(r, n, s, t, a, p, Ftable)
                matches += len(recs)
                points.extend(recs)
    return dict(p=p, r=r % p, n=n % p,
                disc_char=legendre(r*r-4*n, p),
                smooth_base_curves=curves, open_matches=matches,
                points=points, samples=points[:sample_limit])


def scan_fixed_box(primes: Sequence[int], bound: int):
    """Sieve integral B=x^2+r*x+n, |r|,|n|<=bound, disc nonsquare over Q."""
    tables = {p: build_F_hash(p)[0] for p in primes}
    masks = {}
    for p in primes:
        mask = set()
        table = tables[p]
        for r in range(p):
            for n in range(p):
                if fixed_B_census(p, r, n, 0, table)["open_matches"]:
                    mask.add((r, n))
        masks[p] = mask

    survivors = []
    checked = 0
    for r in range(-bound, bound+1):
        for n in range(-bound, bound+1):
            disc = r*r-4*n
            if disc == 0 or (disc > 0 and int(disc**0.5)**2 == disc):
                continue
            checked += 1
            if all((r % p, n % p) in masks[p] for p in primes):
                survivors.append((r, n, disc))
    return dict(primes=list(primes), bound=bound, checked=checked,
                mask_sizes={p: len(masks[p]) for p in primes},
                survivors=survivors)


def print_sample(rec):
    b,w,r,n,s,t,a,tau,q,B,A,S,F = rec
    print(" SAMPLE b,w,r,n,s,t,a=%s tau=%d" %
          ((b,w,r,n,s,t,a),tau))
    print("  B=%s q=%s A=%s" % (B,q,A))
    print("  S=%s F=%s" % (S,F))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", choices=["full","fixed","box"], default="full")
    ap.add_argument("--primes", default="7,11,13")
    ap.add_argument("--discriminant", choices=["nonsplit","split","separable","all"],
                    default="nonsplit")
    ap.add_argument("--r", type=int, default=0)
    ap.add_argument("--n", type=int, default=1)
    ap.add_argument("--bound", type=int, default=20)
    ap.add_argument("--sample-limit", type=int, default=4)
    args = ap.parse_args()
    primes = [int(z) for z in args.primes.split(",") if z]
    if args.mode == "box":
        rep = scan_fixed_box(primes, args.bound)
        print(" ".join(f"{k}={v}" for k,v in rep.items() if k != "survivors"))
        print("survivors=%s" % (rep["survivors"],))
        return
    for p in primes:
        if p in (2,5):
            continue
        if args.mode == "fixed":
            rep = fixed_B_census(p,args.r,args.n,args.sample_limit)
        else:
            rep = full_census(p,args.discriminant,args.sample_limit)
        samples = rep.pop("samples")
        points = rep.pop("points",None)
        print(" ".join(f"{k}={v}" for k,v in rep.items()))
        for rec in samples:
            print_sample(rec)


if __name__ == "__main__":
    main()
