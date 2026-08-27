#!/usr/bin/env python3
"""Lift fresh p=11,13 boundary incidences to explicit open Q_p branches.

Variables are (a,b,c,d,u,t,w,y0,y1,y2,y3).  The eight equations are the
four direct-contact coefficient identities and y_i^2=R_i for the four full
A(2,2,2,8) radicands.  Digit-by-digit multivariate Hensel lifting retains
the square-root variables, so a reported row is stronger than the projected
first-order masks used by the CRT search.
"""

from __future__ import annotations
import argparse,csv,itertools,random
from pathlib import Path

ROOT=Path("results")

def inv(a,m): return pow(a%m,-1,m)

def coeff_equations(x,mod):
    a,b,c,d,u,t,w=x[:7]; aa=[z*z%mod for z in (a,b,c,d)]
    e1=sum(aa)%mod
    e2=sum(aa[i]*aa[j] for i in range(4) for j in range(i+1,4))%mod
    e3=sum(aa[i]*aa[j]*aa[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4))%mod
    e4=aa[0]*aa[1]*aa[2]*aa[3]%mod;i2=inv(2,mod);i4=inv(4,mod)
    s1=(-3*i4*u*u+3*i2*u-3*t*t+2*w+i4)%mod
    s2=(-u**3-6*u*t*t+2*t**3+(1+3*u)*w)%mod
    s3=(-3*u*u*t*t+(1+3*u)*t**3-3*t**4+w*w)%mod
    s4=(-3*u*t**4+2*t**3*w)%mod
    return [(e1-s1)%mod,(e2-s2)%mod,(e3-s3)%mod,(e4-s4)%mod]

def radicands(x,mod):
    a,b,c,d=x[:4]
    return [(a*b*c*d)%mod,
            (a*(a+b)*(a+c)*(a+d))%mod,
            (b*(b+a)*(b+c)*(b+d))%mod,
            (c*(c+a)*(c+b)*(c+d))%mod]

def equations(x,mod):
    return coeff_equations(x,mod)+[(x[7+i]*x[7+i]-r)%mod for i,r in enumerate(radicands(x,mod))]

def jacobian(x,p):
    weights=[-sum(inv(j,p) for j in range(1,5))%p,
             4%p,(-3)%p,(4*inv(3,p))%p,(-inv(4,p))%p]
    J=[[] for _ in range(8)]
    for j in range(11):
        col=[0]*8
        for k,w in enumerate(weights):
            y=list(x);y[j]=(y[j]+k)%p;v=equations(y,p)
            for i in range(8):col[i]=(col[i]+w*v[i])%p
        for i in range(8):J[i].append(col[i])
    return J

def affine_solutions(A,rhs,p):
    """Return one solution and a nullspace basis for A*x=rhs."""
    M=[[(z%p) for z in row]+[rhs[i]%p] for i,row in enumerate(A)]
    piv=[];r=0;n=len(A[0])
    for c in range(n):
        k=next((i for i in range(r,len(M)) if M[i][c]),None)
        if k is None:continue
        M[r],M[k]=M[k],M[r];z=inv(M[r][c],p);M[r]=[(q*z)%p for q in M[r]]
        for i in range(len(M)):
            if i!=r and M[i][c]:
                z=M[i][c];M[i]=[(M[i][j]-z*M[r][j])%p for j in range(n+1)]
        piv.append(c);r+=1
    if any(not any(row[:n]) and row[n] for row in M):return None,[],r
    free=[j for j in range(n) if j not in piv];part=[0]*n
    for i,c in enumerate(piv):part[c]=M[i][n]
    basis=[]
    for f in free:
        v=[0]*n;v[f]=1
        for i,c in enumerate(piv):v[c]=(-M[i][f])%p
        basis.append(v)
    return part,basis,r

def add_solution(part,basis,cs,p):
    return [(part[j]+sum(cs[i]*basis[i][j] for i in range(len(basis))))%p for j in range(len(part))]

def degeneracy_values(x,mod):
    a,b,c,d,u,t,w=x[:7];bs=[a,b,c,d]
    vals={"t":t,"qdisc":u*u-4*t*t,"rootprod":a*b*c*d}
    for i in range(4):
        vals[f"qroot{i}"]=bs[i]**4-u*bs[i]**2+t*t
        vals[f"y{i}"]=x[7+i]
    for i in range(4):
        for j in range(i+1,4): vals[f"diff{i}{j}"]=bs[i]**2-bs[j]**2
    vals["rect01_23"]=a*b-c*d;vals["rect02_13"]=a*c-b*d;vals["rect03_12"]=a*d-b*c
    return {k:v%mod for k,v in vals.items()}

def is_open(x,mod):
    # Smooth/direct-chart openness.  Rectangle equalities are logged separately:
    # they detect a familiar bielliptic locus but are not singularities of the
    # cover/contact equations themselves.
    z=degeneracy_values(x,mod)
    return all(v for k,v in z.items() if not k.startswith("rect"))

def roots_mod(r,p):return [z for z in range(p) if z*z%p==r%p]

def lift_one(base,p,target,rng):
    rads=radicands(base,p);rootsets=[roots_mod(z,p) for z in rads]
    starts=[list(base)+list(ys) for ys in itertools.product(*rootsets)]
    assert starts and all(all(v==0 for v in equations(x,p)) for x in starts)
    jcache={}
    def J(x):
        k=tuple(z%p for z in x)
        if k not in jcache:jcache[k]=jacobian(x,p)
        return jcache[k]
    ranks=sorted(set(affine_solutions(J(x),[0]*8,p)[2] for x in starts))
    states=starts
    for exp in range(1,target):
        pk=p**exp;mod=pk*p;new={}
        for x in states:
            vals=equations(x,mod)
            if any(v%pk for v in vals):continue
            rhs=[-(v//pk)%p for v in vals]
            part,basis,rank=affine_solutions(J(x),rhs,p)
            if part is None:continue
            # Exhaust the small relative tangent space on the first digit;
            # thereafter sample enough continuations for singular branches.
            total=p**len(basis)
            if total<=3000 and len(states)<=32:
                coeffs=itertools.product(range(p),repeat=len(basis))
            else:
                coeffs=([rng.randrange(p) for _ in basis] for __ in range(160))
            for cs in coeffs:
                delta=add_solution(part,basis,cs,p)
                y=tuple((x[j]+pk*delta[j])%mod for j in range(11))
                new[y]=sum(v!=0 for v in degeneracy_values(y,mod).values())
        if not new:return None,ranks
        # On rank-deficient boundary components, many apparently open p^2
        # points die at p^3.  Perform one exact consistency look-ahead before
        # pruning the beam, so those second-order obstructions are respected.
        if exp+1<target:
            nextpk=mod;nextmod=mod*p;liftable={}
            for y,score in new.items():
                vv=equations(y,nextmod)
                if any(v%nextpk for v in vv):continue
                rr=[-(v//nextpk)%p for v in vv]
                part2,basis2,rank2=affine_solutions(J(y),rr,p)
                if part2 is not None:liftable[y]=score
            new=liftable
            if not new:return None,ranks
        opened=[x for x in new if is_open(x,mod)]
        if opened:
            # Once a quantity is nonzero modulo p^k it stays nonzero in every
            # continuation.  Keep a small open beam for robust further lifts.
            # These boundary points can have Jacobian rank seven, hence a
            # genuine second-order obstruction.  Retain enough distinct open
            # p^k states that the next 1/p consistency layer is represented.
            states=opened[:512]
        else:
            states=[x for x,s in sorted(new.items(),key=lambda kv:-kv[1])[:160]]
    ans=next((x for x in states if is_open(x,p**target)),None)
    return ans,ranks

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--exponent",type=int,default=5)
    ap.add_argument("--seed",type=int,default=222245);ap.add_argument("--output")
    ap.add_argument("--log");ap.add_argument("--primes",default="11,13")
    ap.add_argument("--max-rows",type=int,default=0)
    ap.add_argument("--base",help="optional comma-separated a,b,c,d,u,t,w boundary incidence")
    args=ap.parse_args();rng=random.Random(args.seed)
    out=Path(args.output or ROOT/"target_22224_direct_contact_deep13_padic_lift.tsv")
    log=Path(args.log or ROOT/"target_22224_direct_contact_deep13_padic_lift.log")
    found=[];lines=[]
    for p in (int(z) for z in args.primes.split(",") if z):
        inp=ROOT/f"target_22224_direct_contact_deep13_boundary_p{p}.tsv"
        rows=[]
        with inp.open() as f:
            for r in csv.DictReader(f,delimiter="\t"):
                base=tuple(int(r[k])%p for k in ("a","b","c","d","u","t","w"))
                rz=radicands(base,p);ys=[roots_mod(z,p)[0] for z in rz]
                rank=affine_solutions(jacobian(list(base)+ys,p),[0]*8,p)[2]
                rows.append((base,r["flags"],rank,sum(z==0 for z in rz)))
        # Smoothest explicit incidence components first.  In particular this
        # reaches the rank-eight zero-radicand branches at p=13 before the
        # highly singular all-zero rows appearing first in the TSV.
        rows.sort(key=lambda z:(-z[2],z[3]==0,z[3]))
        if args.base:
            wanted=tuple(int(z)%p for z in args.base.split(","));rows=[z for z in rows if z[0]==wanted]
        if args.max_rows:rows=rows[:args.max_rows]
        cert=None
        for base,flags,rank0,nzero in rows:
            x,ranks=lift_one(base,p,args.exponent,rng)
            if x is not None:
                cert=(base,flags,x,ranks);break
        if cert is None:
            lines.append(f"PADIC_LIFT p {p} NO_OPEN_CERT tried {len(rows)}")
            continue
        base,flags,x,ranks=cert;mod=p**args.exponent;eq=equations(x,mod);deg=degeneracy_values(x,mod)
        assert all(v==0 for v in eq) and is_open(x,mod)
        found.append((p,args.exponent,base,flags,x))
        lines.append(f"PADIC_LIFT p {p} exponent {args.exponent} modulus {mod} base {base} flags {flags} jacobian_ranks {ranks}")
        lines.append(f"  vector {x}")
        lines.append(f"  equations {eq}")
        lines.append(f"  open_values {deg}")
    with out.open("w") as f:
        f.write("p\texponent\tbase_a\tbase_b\tbase_c\tbase_d\tbase_u\tbase_t\tbase_w\tflags\ta\tb\tc\td\tu\tt\tw\ty0\ty1\ty2\ty3\n")
        for p,e,base,flags,x in found:
            f.write("\t".join(map(str,(p,e,*base,flags,*x)))+"\n")
    lines.append(f"CERTIFICATES {len(found)} OUTPUT {out}")
    text="\n".join(lines)+"\n";print(text,end="");log.write_text(text)

if __name__=="__main__":main()
