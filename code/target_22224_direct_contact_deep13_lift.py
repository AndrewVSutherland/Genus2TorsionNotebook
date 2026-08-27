#!/usr/bin/env python3
"""First-order lifts of the *fresh* direct-contact boundary incidences.

The input rows come from target_22224_direct_contact_deep13_boundary.m.
For each signed branch incidence modulo p, linearize the four exact
coefficient identities in

    (a,b,c,d,u,t,w).

Whenever a full-cover radicand is zero modulo p, its first-order term must
also vanish: a value of p-adic valuation one cannot be a square.  The output
is the complete projected mask of (u,t,w) modulo p^2 admitted by these
first-order equations.  No superseded Aaux equation is used.
"""

from __future__ import annotations
import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path


def inv(a: int, m: int) -> int:
    return pow(a % m, -1, m)


def coeff_equations(x, mod):
    a, b, c, d, u, t, w = (z % mod for z in x)
    q = [a*a % mod, b*b % mod, c*c % mod, d*d % mod]
    e1 = sum(q) % mod
    e2 = sum(q[i]*q[j] for i in range(4) for j in range(i+1,4)) % mod
    e3 = sum(q[i]*q[j]*q[k] for i in range(4) for j in range(i+1,4)
             for k in range(j+1,4)) % mod
    e4 = q[0]*q[1]*q[2]*q[3] % mod
    i2, i4 = inv(2, mod), inv(4, mod)
    s1 = (-3*i4*u*u + 3*i2*u - 3*t*t + 2*w + i4) % mod
    s2 = (-u**3 - 6*u*t*t + 2*t**3 + (1+3*u)*w) % mod
    s3 = (-3*u*u*t*t + (1+3*u)*t**3 - 3*t**4 + w*w) % mod
    s4 = (-3*u*t**4 + 2*t**3*w) % mod
    return [(e1-s1) % mod, (e2-s2) % mod,
            (e3-s3) % mod, (e4-s4) % mod]


def radicands(x, mod):
    a,b,c,d=x[:4]
    return [(a*b*c*d) % mod,
            (a*(a+b)*(a+c)*(a+d)) % mod,
            (b*(b+a)*(b+c)*(b+d)) % mod,
            (c*(c+a)*(c+b)*(c+d)) % mod]


def gradient(func, x, p):
    # Formal derivative, not the forward difference f(x+1)-f(x).
    # Every equation has degree at most four in each individual variable, so
    # the derivative at zero of the degree-four Lagrange interpolant is exact.
    weights=[-sum(inv(j,p) for j in range(1,5))%p,
             4%p,(-3)%p,(4*inv(3,p))%p,(-inv(4,p))%p]
    ans=[]
    for j in range(7):
        val=0
        for k,w in enumerate(weights):
            y=list(x);y[j]=(y[j]+k)%p
            val=(val+w*func(y,p))%p
        ans.append(val)
    return ans


def consistent(A, rhs, p):
    """Whether A*x=rhs over F_p has a solution."""
    if not A:
        return True
    M=[[(z%p) for z in row]+[rhs[i]%p] for i,row in enumerate(A)]
    r=0;n=4
    for c in range(n):
        piv=next((i for i in range(r,len(M)) if M[i][c]),None)
        if piv is None: continue
        M[r],M[piv]=M[piv],M[r]
        z=inv(M[r][c],p);M[r]=[(v*z)%p for v in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                z=M[i][c];M[i]=[(M[i][j]-z*M[r][j])%p for j in range(n+1)]
        r+=1
    return all(any(row[j] for j in range(n)) or row[n]==0 for row in M)


def main():
    ap=argparse.ArgumentParser();ap.add_argument("--prime",type=int,default=13)
    ap.add_argument("--input");ap.add_argument("--output");ap.add_argument("--log")
    args=ap.parse_args();p=args.prime;p2=p*p
    base="results/target_22224_direct_contact_deep13_boundary"
    inp=Path(args.input or f"{base}_p{p}.tsv")
    out=Path(args.output or f"results/target_22224_direct_contact_deep13_lift_p{p}.tsv")
    log=Path(args.log or f"results/target_22224_direct_contact_deep13_lift_p{p}.log")
    rows=[]
    with inp.open() as f:
        for r in csv.DictReader(f,delimiter="\t"):
            rows.append(tuple(int(r[k])%p for k in ("a","b","c","d","u","t","w")))
    live=set();bybase=defaultdict(set);zero_counts=Counter();inc_live=0;inc_dead=0
    for x in rows:
        # Four coefficient equations plus each zero-radicand equation.
        funcs=[]
        for i in range(4): funcs.append(lambda y,m,i=i: coeff_equations(y,m)[i])
        rz=radicands(x,p)
        for i,z in enumerate(rz):
            if z==0: funcs.append(lambda y,m,i=i: radicands(y,m)[i])
        zero_counts[len(funcs)-4]+=1
        branch_rows=[];param_rows=[];const=[]
        for fn in funcs:
            val=fn(x,p2)%p2
            assert val%p==0
            const.append((val//p)%p)
            g=gradient(fn,x,p)
            branch_rows.append(g[:4]);param_rows.append(g[4:])
        local=set()
        for du in range(p):
            for dt in range(p):
                for dw in range(p):
                    dp=(du,dt,dw)
                    rhs=[(-const[i]-sum(param_rows[i][j]*dp[j] for j in range(3)))%p
                         for i in range(len(funcs))]
                    if consistent(branch_rows,rhs,p):
                        U=x[4]+p*du;T=x[5]+p*dt;W=x[6]+p*dw
                        local.add((U,T,W));live.add((U,T,W));
                        bybase[(x[4],x[5],x[6])].add((U,T,W))
        if local:inc_live+=1
        else:inc_dead+=1
    with out.open("w") as f:
        f.write("u\tt\tw\tu0\tt0\tw0\n")
        for u,t,w in sorted(live):f.write(f"{u}\t{t}\t{w}\t{u%p}\t{t%p}\t{w%p}\n")
    lines=[f"DIRECT_CONTACT_FRESH_LIFT p {p} incidences {len(rows)} incidence_live {inc_live} incidence_dead {inc_dead}",
           f"ZERO_RADICAND_EQUATION_COUNTS {dict(sorted(zero_counts.items()))}",
           f"BASE_TRIPLES {len(bybase)} PROJECTED_P2_TRIPLES {len(live)}",
           "PER_BASE "+str(sorted((k,len(v)) for k,v in bybase.items())),
           "SCOPE complete first-order lift of fresh coefficient identities plus necessary zero-radicand valuation parity",
           f"OUTPUT {out}"]
    text="\n".join(lines)+"\n";print(text,end="");log.write_text(text)

if __name__=="__main__":main()
