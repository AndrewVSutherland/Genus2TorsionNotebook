#!/usr/bin/env python3
"""Memory-bounded root-pair census for the full quadratic-B M(12)+5 cover.

Write the two distinct rational roots of B as

    x = -c, d-c,              x = d*Z-c,
    B = b2*(x+c)*(x+c-d),     d != 0.

After scaling A and q by d^5 and d^2, respectively, the norm identity is

    Atilde^2 - lambda^2*Z^2*(Z-1)^2*G = Q^5,
    G=F(d*Z-c), lambda=b2/d^3.

At Z=0,1 put Q(0)=s^2 and Q(1)=t^2.  The signs of s and t absorb the
four possible signs of A at these roots.  Hermite interpolation fixes four
of the five lower coefficients of the monic quintic Atilde; only its Z^4
coefficient ``a`` remains.  Hence

    S=(Atilde^2-Q^5)/(Z^2*(Z-1)^2)

is an explicit quintic in (s,t,a), and the full norm condition is simply
that S and G are proportional by a nonzero square.

This script hashes the two projective quintics.  The exhaustive full-surface
census therefore costs O(p^3+p^4), not O(p^7), and stores only small tuples.
It also audits the genuinely quadratic fixed-root slices {0,-1}, {0,Lroot},
and {-1,Lroot}.  It uses no computer algebra packages and negligible memory.
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from dataclasses import dataclass
from typing import Iterable, Sequence


def trim(f: list[int]) -> list[int]:
    while len(f) > 1 and f[-1] == 0:
        f.pop()
    return f


def add(f: Sequence[int], g: Sequence[int], p: int) -> list[int]:
    n = max(len(f), len(g))
    return trim([((f[i] if i < len(f) else 0) +
                  (g[i] if i < len(g) else 0)) % p for i in range(n)])


def neg(f: Sequence[int], p: int) -> list[int]:
    return [(-x) % p for x in f]


def sub(f: Sequence[int], g: Sequence[int], p: int) -> list[int]:
    return add(f, neg(g, p), p)


def scale(f: Sequence[int], a: int, p: int) -> list[int]:
    return trim([(a*x) % p for x in f])


def mul(f: Sequence[int], g: Sequence[int], p: int) -> list[int]:
    h = [0]*(len(f)+len(g)-1)
    for i, x in enumerate(f):
        for j, y in enumerate(g):
            h[i+j] = (h[i+j] + x*y) % p
    return trim(h)


def power(f: Sequence[int], n: int, p: int) -> list[int]:
    ans = [1]
    base = list(f)
    while n:
        if n & 1:
            ans = mul(ans, base, p)
        base = mul(base, base, p)
        n >>= 1
    return ans


def derivative(f: Sequence[int], p: int) -> list[int]:
    if len(f) <= 1:
        return [0]
    return trim([(i*f[i]) % p for i in range(1, len(f))])


def divmod_poly(f: Sequence[int], g: Sequence[int], p: int):
    r = trim(list(f))
    g = trim(list(g))
    if g == [0]:
        raise ZeroDivisionError
    if len(r) < len(g):
        return [0], r
    q = [0]*(len(r)-len(g)+1)
    gi = pow(g[-1], -1, p)
    while r != [0] and len(r) >= len(g):
        k = len(r)-len(g)
        a = r[-1]*gi % p
        q[k] = a
        for j, gj in enumerate(g):
            r[j+k] = (r[j+k]-a*gj) % p
        trim(r)
    return trim(q), trim(r)


def gcd_poly(f: Sequence[int], g: Sequence[int], p: int) -> list[int]:
    f, g = trim(list(f)), trim(list(g))
    while g != [0]:
        _, r = divmod_poly(f, g, p)
        f, g = g, r
    if f == [0]:
        return f
    return scale(f, pow(f[-1], -1, p), p)


def squarefree(f: Sequence[int], p: int) -> bool:
    return len(gcd_poly(f, derivative(f, p), p)) == 1


def pad(f: Sequence[int], n: int) -> tuple[int, ...]:
    return tuple(list(f)+[0]*(n-len(f)))


def projective_key(f: Sequence[int], p: int):
    """Return (normalized coefficients, pivot value), low degree first."""
    for x in f:
        if x % p:
            inv = pow(x, -1, p)
            return tuple((inv*y) % p for y in f), x % p
    return None, 0


def is_square_nonzero(x: int, p: int) -> bool:
    return x % p != 0 and pow(x % p, (p-1)//2, p) == 1


def square_root(x: int, p: int) -> int:
    # Census primes are small; deterministic root is useful in samples.
    for y in range(1, p):
        if y*y % p == x % p:
            return y
    raise ValueError("not a nonzero square")


def hermite_S(s: int, t: int, a: int, p: int) -> tuple[list[int], list[int], list[int]]:
    """Return S, Atilde, Q for the signed endpoint roots s,t."""
    inv2 = pow(2, -1, p)
    s2, s3, s5 = s*s % p, pow(s, 3, p), pow(s, 5, p)
    t2, t3, t5 = t*t % p, pow(t, 3, p), pow(t, 5, p)
    e = (t2-s2-1) % p
    q = [s2, e, 1]

    c0 = s5
    c1 = 5*inv2*s3*e % p
    c2 = (a + 2*s5 - 5*s3*t2 + 5*s3 +
          5*inv2*s2*t3 + inv2*t5 - 5*inv2*t3 + 2) % p
    c3 = (-2*a - inv2*s5 + 5*inv2*s3*t2 - 5*inv2*s3 -
          5*inv2*s2*t3 + inv2*t5 + 5*inv2*t3 - 3) % p
    at = [c0, c1, c2, c3, a % p, 1]

    numerator = sub(mul(at, at, p), power(q, 5, p), p)
    quotient, remainder = divmod_poly(numerator, [0, 0, 1, -2 % p, 1], p)
    assert remainder == [0]
    assert len(quotient) <= 6
    return list(pad(quotient, 6)), at, q


def compact_F_transformed(b: int, w: int, c: int, d: int, p: int) -> list[int]:
    """G(Z)=F_{b,w}(dZ-c), coefficients low degree first."""
    x = [(-c) % p, d % p]
    one = [1]
    L = add([b % p], scale(x, (2*b-1) % p, p), p)
    H = add(x, scale(add(one, scale(x, b, p), p), w, p), p)
    term1 = mul(L, mul(H, H, p), p)
    onepx = add(one, x, p)
    wL_minus_x2 = sub(scale(L, w, p), mul(x, x, p), p)
    term2 = scale(mul(mul(onepx, onepx, p), wL_minus_x2, p), 4*b, p)
    return list(pad(mul(L, add(term1, term2, p), p), 6))


def residual_values(vals: Sequence[int], p: int) -> list[int]:
    b, w, c, d, s, t, a = vals
    S, _, _ = hermite_S(s, t, a, p)
    G = compact_F_transformed(b, w, c, d, p)
    # lc(G)=-4*b*(2*b-1)*d^5 is nonzero on the chart.
    return [(S[i]*G[5]-S[5]*G[i]) % p for i in range(5)]


@dataclass
class Dual:
    v: int
    g: list[int]
    p: int

    def __add__(self, other):
        other = as_dual(other, self.p, len(self.g))
        return Dual((self.v+other.v) % self.p,
                    [(x+y) % self.p for x, y in zip(self.g, other.g)], self.p)

    __radd__ = __add__

    def __neg__(self):
        return Dual(-self.v % self.p, [(-x) % self.p for x in self.g], self.p)

    def __sub__(self, other):
        return self + (-as_dual(other, self.p, len(self.g)))

    def __rsub__(self, other):
        return as_dual(other, self.p, len(self.g)) - self

    def __mul__(self, other):
        other = as_dual(other, self.p, len(self.g))
        return Dual(self.v*other.v % self.p,
                    [(self.v*y+other.v*x) % self.p
                     for x, y in zip(self.g, other.g)], self.p)

    __rmul__ = __mul__

    def __pow__(self, n: int):
        if n == 0:
            return as_dual(1, self.p, len(self.g))
        ans = as_dual(1, self.p, len(self.g))
        base = self
        while n:
            if n & 1:
                ans = ans*base
            base = base*base
            n >>= 1
        return ans


def as_dual(x, p: int, n: int) -> Dual:
    if isinstance(x, Dual):
        return x
    return Dual(int(x) % p, [0]*n, p)


def dpoly_add(f, g):
    n = max(len(f), len(g)); p=f[0].p if f else g[0].p; ng=len((f or g)[0].g)
    zero=as_dual(0,p,ng)
    return [(f[i] if i<len(f) else zero)+(g[i] if i<len(g) else zero) for i in range(n)]


def dpoly_scale(f, a):
    return [a*x for x in f]


def dpoly_mul(f, g):
    p=f[0].p; ng=len(f[0].g); h=[as_dual(0,p,ng) for _ in range(len(f)+len(g)-1)]
    for i,x in enumerate(f):
        for j,y in enumerate(g): h[i+j]=h[i+j]+x*y
    return h


def dpoly_pow(f,n):
    p=f[0].p; ng=len(f[0].g); ans=[as_dual(1,p,ng)]; base=f
    while n:
        if n&1: ans=dpoly_mul(ans,base)
        base=dpoly_mul(base,base); n>>=1
    return ans


def dual_residuals(vals: Sequence[int], p: int) -> list[Dual]:
    """Automatic derivatives of the five defining equations."""
    n=7
    X=[]
    for i,v in enumerate(vals):
        g=[0]*n; g[i]=1; X.append(Dual(v%p,g,p))
    b,w,c,d,s,t,a=X; inv2=pow(2,-1,p)
    s2=s*s; s3=s**3; s5=s**5; t2=t*t; t3=t**3; t5=t**5
    e=t2-s2-1
    Q=[s2,e,as_dual(1,p,n)]
    A=[s5,
       5*inv2*s3*e,
       a+2*s5-5*s3*t2+5*s3+5*inv2*s2*t3+inv2*t5-5*inv2*t3+2,
       -2*a-inv2*s5+5*inv2*s3*t2-5*inv2*s3-5*inv2*s2*t3+inv2*t5+5*inv2*t3-3,
       a,as_dual(1,p,n)]
    N=dpoly_add(dpoly_mul(A,A),dpoly_scale(dpoly_pow(Q,5),-1))
    # Divide by monic z^4-2z^3+z^2, retaining Dual coefficients.
    rem=N[:]; quo=[as_dual(0,p,n) for _ in range(6)]
    den=[as_dual(0,p,n),as_dual(0,p,n),as_dual(1,p,n),as_dual(-2,p,n),as_dual(1,p,n)]
    for deg in range(9,3,-1):
        qv=rem[deg]; k=deg-4; quo[k]=qv
        for j in range(5): rem[j+k]=rem[j+k]-qv*den[j]
    assert all(x.v==0 and not any(x.g) for x in rem[:4])

    xx=[-c,d]; one=[as_dual(1,p,n)]
    L=dpoly_add([b],dpoly_scale(xx,2*b-1))
    H=dpoly_add(xx,dpoly_scale(dpoly_add(one,dpoly_scale(xx,b)),w))
    term1=dpoly_mul(L,dpoly_mul(H,H))
    onepx=dpoly_add(one,xx)
    wx=dpoly_add(dpoly_scale(L,w),dpoly_scale(dpoly_mul(xx,xx),-1))
    term2=dpoly_scale(dpoly_mul(dpoly_mul(onepx,onepx),wx),4*b)
    G=dpoly_mul(L,dpoly_add(term1,term2))
    while len(G)<6: G.append(as_dual(0,p,n))
    return [quo[i]*G[5]-quo[5]*G[i] for i in range(5)]


def matrix_rank(rows: list[list[int]], p: int) -> int:
    a=[row[:] for row in rows]; m=len(a); n=len(a[0]) if m else 0; r=0
    for j in range(n):
        piv=next((i for i in range(r,m) if a[i][j]%p),None)
        if piv is None: continue
        a[r],a[piv]=a[piv],a[r]
        inv=pow(a[r][j]%p,-1,p); a[r]=[(inv*x)%p for x in a[r]]
        for i in range(m):
            if i!=r and a[i][j]%p:
                q=a[i][j]%p; a[i]=[(x-q*y)%p for x,y in zip(a[i],a[r])]
        r+=1
        if r==m: break
    return r


def smooth_rank(vals: Sequence[int], p: int) -> int:
    rr=dual_residuals(vals,p)
    assert all(x.v==0 for x in rr)
    return matrix_rank([x.g for x in rr],p)


def build_S_hash(p: int):
    table=defaultdict(list); admissible=0
    for s in range(1,p):
        for t in range(1,p):
            for a in range(p):
                S,At,Q=hermite_S(s,t,a,p)
                if not squarefree(Q,p): continue
                key,pivot=projective_key(S,p)
                if key is None or S[5]==0: continue
                admissible+=1
                table[key].append((s,t,a,pivot,Q,At,S))
    return table,admissible


def base_open(b: int,w: int,p: int) -> bool:
    return b%p and w%p and (b-1)%p and (2*b-1)%p


def check_match(b,w,c,d,record,G,gpivot,p):
    s,t,a,spivot,Q,At,S=record
    tau=S[5]*pow(G[5],-1,p)%p
    if not is_square_nonzero(tau,p): return None
    if len(gcd_poly(Q,G,p)) != 1: return None
    vals=(b,w,c,d,s,t,a)
    rank=smooth_rank(vals,p)
    return vals,tau,square_root(tau,p),rank,Q,At,S,G


def full_census(p: int, sample_limit: int):
    table,nS=build_S_hash(p)
    curves=matches=openpts=smooth=signed=0; samples=[]
    for b in range(p):
      for w in range(p):
        if not base_open(b,w,p): continue
        # Smoothness is invariant under the affine root change; test c=0,d=1.
        F=compact_F_transformed(b,w,0,1,p)
        if not squarefree(F,p): continue
        curves+=1
        for c in range(p):
          for d in range(1,p):
            G=compact_F_transformed(b,w,c,d,p)
            key,gpivot=projective_key(G,p)
            if key not in table: continue
            for rec in table[key]:
                matches+=1
                ans=check_match(b,w,c,d,rec,G,gpivot,p)
                if ans is None: continue
                openpts+=1; signed+=2
                if ans[3]==5: smooth+=1
                if len(samples)<sample_limit: samples.append(ans)
    return dict(p=p,S_admissible=nS,S_keys=len(table),smooth_base_curves=curves,
                projective_matches=matches,open_quotient=openpts,
                signed_points=signed,smooth_open=smooth,samples=samples)


def slice_pairs(kind: str,b: int,p: int) -> Iterable[tuple[int,int]]:
    if kind=="zero_minusone":
        yield 0,(-1)%p                 # roots 0,-1
    elif kind=="zero_one":
        yield 0,1                      # roots 0,1
    elif kind=="zero_L":
        yield 0,(-b*pow((2*b-1)%p,-1,p))%p
    elif kind=="minusone_L":
        # root0=-1 means c=1; d=root1-root0=Lroot+1.
        lr=(-b*pow((2*b-1)%p,-1,p))%p
        yield 1,(lr+1)%p
    elif kind=="L_free":
        c=b*pow((2*b-1)%p,-1,p)%p
        for d in range(1,p):
            yield c,d
    else:
        raise ValueError(kind)


def slice_census(p: int,kind: str,sample_limit:int):
    table,nS=build_S_hash(p)
    curves=matches=openpts=smooth=signed=0;samples=[]
    for b in range(p):
      for w in range(p):
        if not base_open(b,w,p): continue
        F=compact_F_transformed(b,w,0,1,p)
        if not squarefree(F,p): continue
        curves+=1
        for c,d in slice_pairs(kind,b,p):
            if d%p==0: continue
            G=compact_F_transformed(b,w,c,d,p)
            key,gpivot=projective_key(G,p)
            if key not in table: continue
            for rec in table[key]:
                matches+=1
                ans=check_match(b,w,c,d,rec,G,gpivot,p)
                if ans is None: continue
                openpts+=1;signed+=2
                # Rank 5 is smooth on the ambient full surface.  A fixed-root
                # zero-dimensional section is transverse iff its 5x5 Jacobian
                # in (b,w,s,t,a) has rank 5; record both through the sample.
                if ans[3]==5:smooth+=1
                if len(samples)<sample_limit:samples.append(ans)
    return dict(p=p,slice=kind,S_admissible=nS,S_keys=len(table),
                smooth_base_curves=curves,projective_matches=matches,
                open_quotient=openpts,signed_points=signed,
                smooth_in_full_surface=smooth,samples=samples)


def print_report(rep):
    samples=rep.pop("samples")
    print(" ".join(f"{k}={v}" for k,v in rep.items()))
    for ans in samples:
        vals,tau,lam,rank,Q,At,S,G=ans
        print(" SAMPLE b,w,c,d,s,t,a=%s tau=%d lambda=%d jacrank=%d"%
              (vals,tau,lam,rank))
        print("  Q=%s Atilde=%s"%(Q,At))
        print("  S=%s G=%s"%(S,G))


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--primes",default="7,11,13")
    ap.add_argument("--mode",choices=["full","slices","all"],default="all")
    ap.add_argument("--sample-limit",type=int,default=2)
    args=ap.parse_args()
    primes=[int(x) for x in args.primes.split(",") if x]
    for p in primes:
        if p in (2,5): continue
        if args.mode in ("full","all"):
            print_report(full_census(p,args.sample_limit))
        if args.mode in ("slices","all"):
            for kind in ("zero_minusone","zero_one","zero_L","minusone_L","L_free"):
                print_report(slice_census(p,kind,args.sample_limit))


if __name__=="__main__":
    main()
