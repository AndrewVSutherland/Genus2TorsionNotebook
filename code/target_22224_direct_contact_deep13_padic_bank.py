#!/usr/bin/env python3
"""Build a small bank of independent explicit open p^5 incidence lifts."""
from __future__ import annotations
import argparse,csv,random
from pathlib import Path
import target_22224_direct_contact_deep13_padic_lift as H

ROOT=Path("results")
ap=argparse.ArgumentParser();ap.add_argument("--prime",type=int,required=True)
ap.add_argument("--count",type=int,default=6);ap.add_argument("--seed",type=int,default=20260718)
args=ap.parse_args();p=args.prime;rng=random.Random(args.seed+p);candidates=[]
with (ROOT/f"target_22224_direct_contact_deep13_boundary_p{p}.tsv").open() as f:
    for r in csv.DictReader(f,delimiter="\t"):
        b=tuple(int(r[k])%p for k in ("a","b","c","d","u","t","w"));rz=H.radicands(b,p)
        ys=[H.roots_mod(z,p)[0] for z in rz]
        rank=H.affine_solutions(H.jacobian(list(b)+ys,p),[0]*8,p)[2];nzero=sum(z==0 for z in rz)
        # p=11 has no rank-eight incidence; zero-radicand rank-seven branches
        # are the ones observed to pass the higher obstruction layers.
        good=(p==11 and rank==7 and nzero>0) or (p==13 and rank==8)
        if good:candidates.append((nzero,b,r["flags"],rank))
candidates.sort(key=lambda z:(z[0],z[1]))
found=[]
for nzero,b,flags,rank in candidates:
    x,ranks=H.lift_one(b,p,5,rng)
    print("BANK_TRY",p,b,"nzero",nzero,"rank",rank,"success",x is not None,flush=True)
    if x is None:continue
    assert all(v==0 for v in H.equations(x,p**5)) and H.is_open(x,p**5)
    found.append((b,flags,x,rank,nzero))
    if len(found)>=args.count:break
out=ROOT/f"target_22224_direct_contact_deep13_padic_bank_p{p}.tsv"
with out.open("w") as f:
    f.write("p\texponent\tbase_a\tbase_b\tbase_c\tbase_d\tbase_u\tbase_t\tbase_w\tflags\trank\tnzero\ta\tb\tc\td\tu\tt\tw\ty0\ty1\ty2\ty3\n")
    for b,flags,x,rank,nzero in found:
        f.write("\t".join(map(str,(p,5,*b,flags,rank,nzero,*x)))+"\n")
text=f"PADIC_BANK p {p} requested {args.count} found {len(found)} candidates {len(candidates)} output {out}\n"
print(text,end="");(ROOT/f"target_22224_direct_contact_deep13_padic_bank_p{p}.log").write_text(text)
