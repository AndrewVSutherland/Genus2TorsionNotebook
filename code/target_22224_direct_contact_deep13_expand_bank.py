#!/usr/bin/env python3
"""Union independently certified p^5 parameter/incidence banks."""
import csv
from pathlib import Path
ROOT=Path("results")
lines=[]
for p in (11,13):
    sources=[ROOT/f"target_22224_direct_contact_deep13_padic_tangent_p{p}.tsv"]
    if p==11:sources.append(ROOT/"target_22224_direct_contact_nearmiss_fixedbranch_p11.tsv")
    rows={}
    for src in sources:
        with src.open() as f:
            for r in csv.DictReader(f,delimiter="\t"):
                q=tuple(int(r[k]) for k in ("u","t","w"));rows.setdefault(q,src.stem)
    out=ROOT/f"target_22224_direct_contact_deep13_padic_expanded_p{p}.tsv"
    with out.open("w") as f:
        f.write("u\tt\tw\tsource\n")
        for q in sorted(rows):f.write("\t".join(map(str,q))+"\t"+rows[q]+"\n")
    lines.append(f"PADIC_EXPANDED p {p} sources {[str(z) for z in sources]} projected {len(rows)} output {out}")
text="\n".join(lines)+"\n";print(text,end="")
(ROOT/"target_22224_direct_contact_deep13_padic_expanded.log").write_text(text)
