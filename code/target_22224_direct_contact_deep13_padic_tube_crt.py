#!/usr/bin/env python3
"""Exhaust the 192 generic-open CRT combinations over certified p^5 tubes."""
from __future__ import annotations
import argparse,csv,itertools,math
from pathlib import Path

ROOT=Path("results")

def triples(path,offset=0):
    with Path(path).open() as f:
        r=csv.reader(f,delimiter="\t");next(r)
        return sorted(set(tuple(int(z[offset+j]) for j in range(3)) for z in r if z))

def crt(a,m,b,n):return (a+m*((b-a)*pow(m,-1,n)%n))%(m*n)

def ratrec(a,m):
    B=math.isqrt(m//2);r0,r1=m,a%m;t0,t1=0,1
    while r1>B:
        q=r0//r1;r0,r1=r1,r0-q*r1;t0,t1=t1,t0-q*t1
    if not t1:return None
    n,d=r1,t1
    if d<0:n,d=-n,-d
    g=math.gcd(n,d);n//=g;d//=g
    if abs(n)>B or d>B or math.gcd(d,m)!=1 or (n-a*d)%m:return None
    return n,d

def allowed(p):
    return set(triples(ROOT/f"target_22224_direct_contact_deep13_boundary_p{p}.tsv",1))|set(triples(ROOT/f"target_22224_direct_contact_deep13_p{p}.tsv",1))

def filter_at(q,p,A):
    if any(d%p==0 for n,d in q):return True,"nonunit"
    z=tuple((n*pow(d,-1,p))%p for n,d in q)
    return z in A,"mask"

cert={}
with (ROOT/"target_22224_direct_contact_deep13_padic_lift.tsv").open() as f:
    for r in csv.DictReader(f,delimiter="\t"):
        p=int(r["p"]);e=int(r["exponent"]);cert[p]=(p**e,(int(r["u"]),int(r["t"]),int(r["w"])))
assert 11 in cert and 13 in cert

ap=argparse.ArgumentParser();ap.add_argument("--single-cert",action="store_true");args=ap.parse_args()
if args.single_cert:
    masks=[(cert[11][0],[cert[11][1]]),(cert[13][0],[cert[13][1]])]
    stem="target_22224_direct_contact_deep13_padic_tube_crt"
else:
    masks=[(cert[11][0],triples(ROOT/"target_22224_direct_contact_deep13_padic_tangent_p11.tsv")),
           (cert[13][0],triples(ROOT/"target_22224_direct_contact_deep13_padic_tangent_p13.tsv"))]
    stem="target_22224_direct_contact_deep13_padic_tangent_crt"
for p in (17,19,23,29,31):masks.append((p,triples(ROOT/f"target_22224_direct_contact_deep13_p{p}.tsv",1)))
A={p:allowed(p) for p in (37,41,43)}
M=math.prod(m for m,z in masks);counts={"comb":0,"ratrec":0,"p37":0,"p41":0,"p43":0};rows=[]
for zz in itertools.product(*(z for m,z in masks)):
    counts["comb"]+=1;q=[]
    for j in range(3):
        a=zz[0][j];m=masks[0][0]
        for i in range(1,len(masks)):
            a=crt(a,m,zz[i][j],masks[i][0]);m*=masks[i][0]
        z=ratrec(a,m)
        if z is None:break
        q.append(z)
    if len(q)!=3 or q[1][0]==0:continue
    counts["ratrec"]+=1;ok=True
    for p in (37,41,43):
        yes,why=filter_at(q,p,A[p])
        if not yes:ok=False;break
        counts[f"p{p}"]+=1
    if ok:rows.append(q)

out=ROOT/f"{stem}.tsv"
with out.open("w") as f:
    f.write("u_num\tu_den\tt_num\tt_den\tw_num\tw_den\n")
    for q in rows:f.write("\t".join(str(v) for z in q for v in z)+"\n")
text=(f"DIRECT_CONTACT_PADIC_TUBE_CRT modulus {M} bound {math.isqrt(M//2)} "
      f"mask_sizes {[len(z) for m,z in masks]} counts {counts} survivors {len(rows)} output {out}\n")
print(text,end="");(ROOT/f"{stem}.log").write_text(text)
