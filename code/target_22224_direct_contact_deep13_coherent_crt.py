#!/usr/bin/env python3
"""Reconstruct coherent branch/contact incidences, not only (u,t,w)."""
from __future__ import annotations
import argparse,csv,math
from fractions import Fraction as F
from pathlib import Path

ROOT=Path("results");M11=11**5;M13=13**5
ap=argparse.ArgumentParser();ap.add_argument("--include",default="");args=ap.parse_args()

def read_rows(paths):
    ans={}
    for path in paths:
        with Path(path).open() as f:
            for r in csv.DictReader(f,delimiter="\t"):
                if not r:continue
                z=tuple(int(r[k]) for k in ("a","b","c","d","u","t","w"));ans[z]=z
    return sorted(ans)
def triples(path,offset=0):
    with Path(path).open() as f:
        z=csv.reader(f,delimiter="\t");next(z)
        return set(tuple(int(r[offset+j]) for j in range(3)) for r in z if r)
def allowed(p):
    return triples(ROOT/f"target_22224_direct_contact_deep13_boundary_p{p}.tsv",1)|triples(ROOT/f"target_22224_direct_contact_deep13_p{p}.tsv",1)
def crt(a,m,b,n):return (a+m*((b-a)*pow(m,-1,n)%n))%(m*n)
def ratrec(a):
    r0,r1=M,a%M;t0,t1=0,1
    while r1>B:
        q=r0//r1;r0,r1=r1,r0-q*r1;t0,t1=t1,t0-q*t1
    if not t1:return None
    n,d=r1,t1
    if d<0:n,d=-n,-d
    g=math.gcd(n,d);n//=g;d//=g
    if abs(n)>B or d>B or math.gcd(d,M)!=1 or (n-a*d)%M:return None
    return F(n,d)
def pass_at(q,p,A):
    pars=q[4:]
    if any(x.denominator%p==0 for x in pars):return True
    z=tuple((x.numerator*pow(x.denominator,-1,p))%p for x in pars)
    return z in A
def coeffs(q):
    a,b,c,d,u,t,w=q;aa=[z*z for z in (a,b,c,d)]
    e1=sum(aa);e2=sum(aa[i]*aa[j] for i in range(4) for j in range(i+1,4))
    e3=sum(aa[i]*aa[j]*aa[k] for i in range(4) for j in range(i+1,4) for k in range(j+1,4));e4=math.prod(aa)
    s1=-F(3,4)*u*u+F(3,2)*u-3*t*t+2*w+F(1,4)
    s2=-u**3-6*u*t*t+2*t**3+(1+3*u)*w
    s3=-3*u*u*t*t+(1+3*u)*t**3-3*t**4+w*w
    s4=-3*u*t**4+2*t**3*w
    return (e1-s1,e2-s2,e3-s3,e4-s4)
def squareq(x):
    return x>=0 and math.isqrt(x.numerator)**2==x.numerator and math.isqrt(x.denominator)**2==x.denominator
def rads(q):
    a,b,c,d=q[:4]
    return (a*b*c*d,a*(a+b)*(a+c)*(a+d),b*(b+a)*(b+c)*(b+d),c*(c+a)*(c+b)*(c+d))

p11=read_rows([ROOT/"target_22224_direct_contact_deep13_padic_tangent_p11.tsv",ROOT/"target_22224_direct_contact_nearmiss_fixedbranch_p11.tsv"])
p13=read_rows([ROOT/"target_22224_direct_contact_deep13_padic_tangent_p13.tsv"])
included=[int(z) for z in args.include.split(",") if z];good=[(1,(0,0,0,0,0,0,0))]
for p in included:
    rows=[]
    with (ROOT/f"target_22224_direct_contact_deep13_p{p}.tsv").open() as f:
        for r in csv.DictReader(f,delimiter="\t"):
            rows.append(tuple(int(r[k]) for k in ("a","b","c","d","u","t","w")))
    rows=sorted(set(rows));new=[]
    for m,x in good:
        for y in rows:new.append((m*p,tuple(crt(x[j],m,y[j],p) for j in range(7))))
    good=new
Gmod=good[0][0];M=M11*M13*Gmod;B=math.isqrt(M//2)
AA={p:allowed(p) for p in (17,19,23,29,31,37,41,43)}
counts={"pairs":0,"ratrec7":0,"exact_contact":0,"full_cover":0,"finite_masks":0};hits=[]
for x in p11:
    for y in p13:
        deep=tuple(crt(a,M11,b,M13) for a,b in zip(x,y));MD=M11*M13
        for gm,g in good:
            counts["pairs"]+=1;q=[]
            for a,b in zip(deep,g):
                z=ratrec(crt(a,MD,b,Gmod))
                if z is None:break
                q.append(z)
            if len(q)!=7 or q[5]==0:continue
            counts["ratrec7"]+=1
            if any(coeffs(q)):continue
            counts["exact_contact"]+=1
            if not all(squareq(z) and z for z in rads(q)):continue
            counts["full_cover"]+=1
            if not all(pass_at(q,p,AA[p]) for p in AA):continue
            counts["finite_masks"]+=1;hits.append(q)

tag="_"+"p"+"p".join(map(str,included)) if included else ""
out=ROOT/f"target_22224_direct_contact_deep13_coherent_crt{tag}.tsv"
names=("a","b","c","d","u","t","w")
with out.open("w") as f:
    f.write("\t".join(v+s for v in names for s in ("_num","_den"))+"\n")
    for q in hits:f.write("\t".join(str(v) for z in q for v in (z.numerator,z.denominator))+"\n")
text=f"DIRECT_CONTACT_COHERENT_CRT p11 {len(p11)} p13 {len(p13)} included {included} good_states {len(good)} modulus {M} bound {B} counts {counts} output {out}\n"
print(text,end="");(ROOT/f"target_22224_direct_contact_deep13_coherent_crt{tag}.log").write_text(text)
