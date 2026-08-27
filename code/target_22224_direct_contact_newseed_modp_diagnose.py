#!/usr/bin/env python3
"""Exact mod-11/13 contact-incidence diagnosis of the newest A2228 seeds."""
import itertools,math
from pathlib import Path
import target_22224_direct_contact_deep13_padic_lift as H
ROOT=Path("results")
seeds={"F4":(-126,28,49,-50),"F4_partner":(-126,28,49,-578),
       "F5":(-112,14,49,-338),"F7":(-18,1,16,-1682),
       "B2000":(144,697,-722,-1394)}
lines=[]
for name,B in seeds.items():
    rr=H.radicands(B,10**40);roots=[]
    for z in rr:
        y=math.isqrt(z);assert y*y==z;roots.append(y)
    lines.append(f"SEED {name} branches {B} radicand_roots {tuple(roots)}")
    for p in (11,13):
        sol=[]
        for lam in range(1,p):
            bb=tuple(lam*z for z in B)
            for u,t,w in itertools.product(range(p),repeat=3):
                if not any(H.coeff_equations((*bb,u,t,w),p)):sol.append((lam,u,t,w))
        signed=tuple(z%p for z in B);squares=tuple(z*z%p for z in signed)
        smooth=all(signed) and len(set(squares))==4
        lines.append(f"  p {p} signed {signed} squares {squares} smooth {smooth} scaled_contact_solutions {sol}")
text="\n".join(lines)+"\n";print(text,end="")
(ROOT/"target_22224_direct_contact_newseed_modp_diagnose.log").write_text(text)
