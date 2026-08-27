#!/usr/bin/env python3
"""Lift the corrected-mask near miss [50,528,-726,-891] along its fixed branch line."""
from __future__ import annotations
import itertools,math
from pathlib import Path
import target_22224_direct_contact_deep13_padic_lift as H

ROOT=Path("results");B=(50,528,-726,-891)
R=H.radicands(B,10**30);Y=tuple(math.isqrt(z) for z in R);assert all(y*y==r for y,r in zip(Y,R))

def fixed11():
    p=11;states=[]
    for u,t,w in itertools.product(range(p),repeat=3):
        if not any(H.coeff_equations((*B,u,t,w),p)):states.append((u,t,w))
    counts=[len(states)]
    for e in range(1,5):
        pk=p**e;mod=pk*p;new=[]
        for q in states:
            for d in itertools.product(range(p),repeat=3):
                z=tuple(q[j]+pk*d[j] for j in range(3))
                if not any(H.coeff_equations((*B,*z),mod)):new.append(z)
        states=new;counts.append(len(states))
    mod=p**5;bb=tuple(z%mod for z in B);yy=tuple(z%mod for z in Y);openstates=[]
    for u,t,w in states:
        x=(*bb,u,t,w,*yy)
        if H.is_open(x,mod):openstates.append((u,t,w,x))
    return counts,openstates

def weights(d,p):
    return [-sum(pow(j,-1,p) for j in range(1,d+1))%p]+[((-1)**(j+1)*math.comb(d,j)*pow(j,-1,p))%p for j in range(1,d+1)]

def fixed13():
    p=13
    def fun(z,m):
        lam,u,t,w=z;bb=tuple(lam*q for q in B)
        return H.coeff_equations((*bb,u,t,w),m)
    def J(z):
        A=[[] for _ in range(4)]
        for var in range(4):
            col=[0]*4
            for k,ww in enumerate(weights(8,p)):
                y=list(z);y[var]=(y[var]+k)%p;v=fun(y,p)
                for i in range(4):col[i]=(col[i]+ww*v[i])%p
            for i in range(4):A[i].append(col[i])
        return A
    states=[]
    for lam in range(1,p):
        for u,t,w in itertools.product(range(p),repeat=3):
            if not any(fun((lam,u,t,w),p)):states.append((lam,u,t,w))
    counts=[len(states)];ranks=sorted(set(H.affine_solutions(J(z),[0]*4,p)[2] for z in states))
    for e in range(1,5):
        pk=p**e;mod=pk*p;new=[]
        for z in states:
            vv=fun(z,mod);rhs=[-(x//pk)%p for x in vv]
            part,basis,rank=H.affine_solutions(J(z),rhs,p)
            if part is None:continue
            for cs in itertools.product(range(p),repeat=len(basis)):
                de=H.add_solution(part,basis,cs,p);y=tuple((z[j]+pk*de[j])%mod for j in range(4))
                assert not any(fun(y,mod));new.append(y)
        states=new;counts.append(len(states))
    return counts,ranks,states

c11,open11=fixed11();c13,ranks13,end13=fixed13()
out=ROOT/"target_22224_direct_contact_nearmiss_fixedbranch_p11.tsv"
with out.open("w") as f:
    f.write("u\tt\tw\tcenter\ta\tb\tc\td\ty0\ty1\ty2\ty3\n")
    for i,(u,t,w,x) in enumerate(open11):f.write("\t".join(map(str,(u,t,w,i,*x[:4],*x[7:])))+"\n")
text=(f"NEARMISS_FIXED_BRANCH branches {B} square_roots {Y}\n"
      f"P11_COUNTS_BY_EXPONENT {c11} open_p5 {len(open11)} output {out}\n"
      f"P13_SCALE_CONTACT_RANKS {ranks13} COUNTS_BY_EXPONENT {c13} terminal {len(end13)}\n"
      "SCOPE fixed exact projective branch ratios; transverse branch deformations are not excluded\n")
print(text,end="");(ROOT/"target_22224_direct_contact_nearmiss_fixedbranch_lift.log").write_text(text)
