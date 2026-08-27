#!/usr/bin/env python3
"""Transverse p13^5 lifts of every fixed-near-miss p13^4 state."""
from __future__ import annotations
import itertools,math,random
from pathlib import Path
import target_22224_direct_contact_deep13_padic_lift as H

ROOT=Path("results");p=13;B=(50,528,-726,-891)
RR=H.radicands(B,10**30);Y=tuple(math.isqrt(z) for z in RR);assert all(y*y==r for y,r in zip(Y,RR))

def weights(d):return [-sum(pow(j,-1,p) for j in range(1,d+1))%p]+[((-1)**(j+1)*math.comb(d,j)*pow(j,-1,p))%p for j in range(1,d+1)]
def fixed_fun(z,m):
    lam,u,t,w=z;bb=tuple(lam*q for q in B)
    return H.coeff_equations((*bb,u,t,w),m)
def fixed_J(z):
    A=[[] for _ in range(4)]
    for var in range(4):
        col=[0]*4
        for k,ww in enumerate(weights(8)):
            y=list(z);y[var]=(y[var]+k)%p;v=fixed_fun(y,p)
            for i in range(4):col[i]=(col[i]+ww*v[i])%p
        for i in range(4):A[i].append(col[i])
    return A
def fixed_states_p4():
    states=[(6,2,1,10),(7,2,1,10)];counts=[len(states)]
    history=[states]
    for e in range(1,4):
        pk=p**e;mod=pk*p;new=[]
        for z in states:
            vv=fixed_fun(z,mod);rhs=[-(x//pk)%p for x in vv]
            part,basis,rank=H.affine_solutions(fixed_J(z),rhs,p)
            if part is None:continue
            for cs in itertools.product(range(p),repeat=len(basis)):
                de=H.add_solution(part,basis,cs,p);y=tuple((z[j]+pk*de[j])%mod for j in range(4))
                assert not any(fixed_fun(y,mod));new.append(y)
        states=new;counts.append(len(states));history.append(states)
    return counts,states,history
def rank3(vs):
    A=[list(v) for v in vs];r=0
    for c in range(3):
        k=next((i for i in range(r,len(A)) if A[i][c]%p),None)
        if k is None:continue
        A[r],A[k]=A[k],A[r];q=pow(A[r][c]%p,-1,p);A[r]=[(z*q)%p for z in A[r]]
        for i in range(len(A)):
            if i!=r and A[i][c]%p:
                q=A[i][c]%p;A[i]=[(A[i][j]-q*A[r][j])%p for j in range(3)]
        r+=1
    return r
def projected_basis(fullbasis):
    chosen=[];projs=[]
    for v in fullbasis:
        q=tuple(v[j] for j in (4,5,6))
        if rank3(projs+[q])>rank3(projs):chosen.append(v);projs.append(q)
    return chosen

counts,states,history=fixed_states_p4();pk=p**4;mod=p**5;params={};live=0;openstates=0;proj_dims={}
for center,z in enumerate(states):
    lam,u,t,w=z;x=tuple(lam*q%pk for q in B)+(u,t,w)+tuple(lam*lam*q%pk for q in Y)
    assert not any(H.equations(x,pk))
    vv=H.equations(x,mod);rhs=[-(q//pk)%p for q in vv]
    part,basis,rank=H.affine_solutions(H.jacobian(x,p),rhs,p)
    if part is None:continue
    live+=1;pb=projected_basis(basis);proj_dims[len(pb)]=proj_dims.get(len(pb),0)+1
    for cs in itertools.product(range(p),repeat=len(pb)):
        de=H.add_solution(part,pb,cs,p);y=tuple((x[j]+pk*de[j])%mod for j in range(11))
        assert not any(H.equations(y,mod))
        if not H.is_open(y,mod):continue
        openstates+=1;params.setdefault((y[4],y[5],y[6]),(center,y))

# The fixed p^4 tube is obstructed.  Move transversely one digit earlier:
# sample the five-dimensional full-incidence p^3 -> p^4 tangent space and
# retain every sample whose next p^5 linear obstruction is soluble.
early_trials=0;early_live=0;rng=random.Random(50528726);p3states=history[2];pk3=p**3;mod4=p**4
for center,z in enumerate(p3states):
    lam,u,t,w=z;x=tuple(lam*q%pk3 for q in B)+(u,t,w)+tuple(lam*lam*q%pk3 for q in Y)
    assert not any(H.equations(x,pk3));vv=H.equations(x,mod4);rhs=[-(q//pk3)%p for q in vv]
    part,basis,rank=H.affine_solutions(H.jacobian(x,p),rhs,p)
    if part is None:continue
    for trial in range(600):
        early_trials+=1;cs=[rng.randrange(p) for _ in basis];de=H.add_solution(part,basis,cs,p)
        y=tuple((x[j]+pk3*de[j])%mod4 for j in range(11));assert not any(H.equations(y,mod4))
        vv=H.equations(y,mod);rhs=[-(q//mod4)%p for q in vv]
        part2,basis2,rank2=H.affine_solutions(H.jacobian(y,p),rhs,p)
        if part2 is None:continue
        early_live+=1;z5=tuple((y[j]+mod4*part2[j])%mod for j in range(11))
        assert not any(H.equations(z5,mod))
        if H.is_open(z5,mod):params.setdefault((z5[4],z5[5],z5[6]),(len(states)+center,z5))

out=ROOT/"target_22224_direct_contact_nearmiss_transverse_p13.tsv"
with out.open("w") as f:
    f.write("u\tt\tw\tcenter\ta\tb\tc\td\ty0\ty1\ty2\ty3\n")
    for q in sorted(params):
        center,y=params[q];f.write("\t".join(map(str,(*q,center,*y[:4],*y[7:])))+"\n")
text=(f"NEARMISS_TRANSVERSE_P13 fixed_counts_p1_to_p4 {counts} fixed_p4 {len(states)} "
      f"full_lift_live {live} projected_dimensions {proj_dims} open_incidence_representatives {openstates} "
      f"early_transverse_trials {early_trials} early_p5_live {early_live} "
      f"projected_params {len(params)} output {out}\n")
print(text,end="");(ROOT/"target_22224_direct_contact_nearmiss_transverse_p13.log").write_text(text)
