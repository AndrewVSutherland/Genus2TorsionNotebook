#!/usr/bin/env python3
"""Resumable CRT/reconstruction sampler for the fresh direct-contact masks."""
from __future__ import annotations
import argparse,csv,math,random
from pathlib import Path

ROOT=Path("results")

def read_triples(path,cols=("u","t","w")):
    out=set()
    with Path(path).open() as f:
        for r in csv.DictReader(f,delimiter="\t"):
            out.add(tuple(int(r[c]) for c in cols))
    return sorted(out)

def allowed_prime(p):
    b=read_triples(ROOT/f"target_22224_direct_contact_deep13_boundary_p{p}.tsv")
    t=read_triples(ROOT/f"target_22224_direct_contact_deep13_p{p}.tsv")
    return sorted(set(b)|set(t)),len(set(b)),len(set(t))

def crt(a,m,b,n):
    return (a+m*((b-a)*pow(m,-1,n)%n))%(m*n)

def ratrec(a,m):
    """Classical symmetric rational reconstruction, or None."""
    bound=math.isqrt(m//2)
    r0,r1=m,a%m;t0,t1=0,1
    while r1>bound:
        q=r0//r1;r0,r1=r1,r0-q*r1;t0,t1=t1,t0-q*t1
    if t1==0:return None
    n,d=r1,t1
    if d<0:n,d=-n,-d
    g=math.gcd(n,d);n//=g;d//=g
    if abs(n)>bound or d>bound or math.gcd(d,m)!=1:return None
    if (n-a*d)%m:return None
    return n,d

def main():
    ap=argparse.ArgumentParser();ap.add_argument("--trials",type=int,default=1_000_000)
    ap.add_argument("--seed",type=int,default=22224);ap.add_argument("--output")
    ap.add_argument("--log")
    ap.add_argument("--anchor-record",action="store_true",
                    help="fix both depth-two masks to the record-centred class (0,0,0)")
    ap.add_argument("--extra-primes",default="",
                    help="comma-separated additional finite-mask primes, e.g. 37,41,43")
    args=ap.parse_args();rng=random.Random(args.seed)
    masks=[]
    m11=read_triples(ROOT/"target_22224_direct_contact_deep13_lift_p11.tsv")
    m13=read_triples(ROOT/"target_22224_direct_contact_deep13_lift_p13.tsv")
    if args.anchor_record:
        assert (0,0,0) in m11 and (0,0,0) in m13
        m11=[(0,0,0)];m13=[(0,0,0)]
    masks.append((121,m11,"p11^2"));masks.append((169,m13,"p13^2"))
    details=[]
    extra=tuple(int(z) for z in args.extra_primes.split(",") if z.strip())
    primes=tuple(dict.fromkeys((17,19,23,29,31)+extra))
    for p in primes:
        z,nb,nt=allowed_prime(p);masks.append((p,z,f"p{p}"));details.append((p,len(z),nb,nt))
    modulus=math.prod(m for m,_,_ in masks);seen=set();rows=[];reconstructed=0
    for _ in range(args.trials):
        chosen=[rng.choice(z) for _,z,_ in masks]
        vals=[]
        for j in range(3):
            a=chosen[0][j];m=masks[0][0]
            for i in range(1,len(masks)):
                a=crt(a,m,chosen[i][j],masks[i][0]);m*=masks[i][0]
            rr=ratrec(a,m)
            if rr is None:break
            vals.append(rr)
        if len(vals)!=3 or vals[1][0]==0:continue
        reconstructed+=1;key=tuple(vals)
        if key in seen:continue
        seen.add(key);rows.append(vals)
    out=Path(args.output or ROOT/f"target_22224_direct_contact_deep13_crt_seed{args.seed}.tsv")
    log=Path(args.log or ROOT/f"target_22224_direct_contact_deep13_crt_seed{args.seed}.log")
    with out.open("w") as f:
        f.write("u_num\tu_den\tt_num\tt_den\tw_num\tw_den\n")
        for z in rows:f.write("\t".join(str(v) for q in z for v in q)+"\n")
    lines=[f"DIRECT_CONTACT_CRT seed {args.seed} trials {args.trials} anchor_record {args.anchor_record} modulus {modulus} reconstruction_bound {math.isqrt(modulus//2)}",
           f"MASK_SIZES p11sq {len(m11)} p13sq {len(m13)} good {details}",
           f"RECONSTRUCTED {reconstructed} UNIQUE {len(rows)}",f"OUTPUT {out}"]
    text="\n".join(lines)+"\n";print(text,end="");log.write_text(text)
if __name__=="__main__":main()
