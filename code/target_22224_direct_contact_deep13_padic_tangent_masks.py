#!/usr/bin/env python3
"""Enumerate the final-digit tangent slices around the certified p^5 branches."""
from __future__ import annotations
import csv,itertools
from pathlib import Path
import target_22224_direct_contact_deep13_padic_lift as H

ROOT=Path("results")
rows={11:[],13:[]}
sources=[ROOT/"target_22224_direct_contact_deep13_padic_lift.tsv"]
sources+=sorted(ROOT.glob("target_22224_direct_contact_deep13_padic_sample_p*_s*.tsv"))
for src in sources:
    with src.open() as f:
        for r in csv.DictReader(f,delimiter="\t"):
            p=int(r["p"]);rows[p].append((src.stem,r))

lines=[]
for p in (11,13):
    params={};states=0;openstates=0;offrect=0;ranks=set();nullities=set()
    for center,(src,r) in enumerate(rows[p]):
        e=int(r["exponent"]);mod=p**e;step=p**(e-1)
        x=tuple(int(r[k]) for k in ("a","b","c","d","u","t","w","y0","y1","y2","y3"))
        part,basis,rank=H.affine_solutions(H.jacobian(x,p),[0]*8,p)
        assert part is not None and all(z==0 for z in part);ranks.add(rank);nullities.add(len(basis))
        for cs in itertools.product(range(p),repeat=len(basis)):
            delta=H.add_solution(part,basis,cs,p)
            y=tuple((x[j]+step*delta[j])%mod for j in range(11));states+=1
            assert all(v==0 for v in H.equations(y,mod))
            if not H.is_open(y,mod):continue
            openstates+=1;d=H.degeneracy_values(y,mod)
            if all(d[k] for k in d if k.startswith("rect")):offrect+=1
            params.setdefault((y[4],y[5],y[6]),(center,y))
    out=ROOT/f"target_22224_direct_contact_deep13_padic_tangent_p{p}.tsv"
    with out.open("w") as f:
        f.write("u\tt\tw\tcenter\ta\tb\tc\td\ty0\ty1\ty2\ty3\n")
        for q in sorted(params):
            center,y=params[q];f.write("\t".join(map(str,(*q,center,*y[:4],*y[7:])))+"\n")
    lines.append(f"PADIC_TANGENT p {p} centers {len(rows[p])} exponent {e} ranks {sorted(ranks)} nullities {sorted(nullities)} states {states} open {openstates} off_rectangle {offrect} projected_params {len(params)} output {out}")
text="\n".join(lines)+"\n";print(text,end="")
(ROOT/"target_22224_direct_contact_deep13_padic_tangent.log").write_text(text)
