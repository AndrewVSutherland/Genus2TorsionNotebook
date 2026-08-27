#!/usr/bin/env python3
"""Project direct-contact boundary incidences to repeated-fiber T classes."""
from __future__ import annotations
import csv
from pathlib import Path

ROOT=Path("results")
FIBERS={
 "P8":((528,-726,-891),50),
 "F1":((-1470,-630,336),25),
 "F2":((-720,20,300),-363),
 "F3":((-612,34,289),-338),
 "F4":((-126,28,49),-50),
 "F5":((-112,14,49),-50),
 "F6":((-50,30,45),-48),
 "F7":((-18,1,16),-50),
 # New square-ratio fibres from the live B=5000 stream.  Denominators
 # have been cleared by one common projective scaling in each row.
 "N1":((-15,50,90),-48),
 "N2":((-2450,-98,3920),4205),
 "N3":((-529,-92,1150),578),
 "N4":((-468,-39,676),768),
 "N5":((-961,-62,1984),1250),
 "N6":((-507,-338,3042),588),
 "N7":((-75,80,120),-72),
 "N8":((-1372,-1078,3430),1375),
 "N9":((-2300,2350,2450),-1081),
 "N10":((-2205,-588,2450),3630),
 "N11":((-960,1800,2535),-1352),
 "N12":((-2156,3850,4235),-3610),
 "O1":((-27,108,162),-98),
 "O2":((-98,588,882),-507),
 "O3":((-147,-98,588),162),
 "O4":((-72,80,120),-75),
 "O5":((-45,-30,48),50),
 "O6":((-405,-27,675),605),
 "O7":((-72,-40,240),75),
 "O8":((-550,1375,1750),-1372),
 "O9":((9,495,891),5),
 "M1":((-48,-15,90),50),
 "M2":((-32,800,1280),-605),
 "M3":((-325,-13,845),405),
 "M4":((-2401,-98,4900),3362),
 "M5":((99,244,4026),6),
 # New live B=10000 fibres at a approximately 779.
 "R1":((-338,676,8450),-529),
 "R2":((-484,-176,539),676),
 "R3":((-289,-272,850),338),
 "R4":((-175,400,1400),-392),
 "R5":((-529,-414,1472),578),
 "R6":((-363,605,9075),-405),
 "R7":((-2209,-188,4606),2738),
 # Final three fibres first visible after a=779 in the completed stream.
 "S1":((-980,2205,7350),-2166),
 # S2 is the same projective family as N11, retained as an explicit alias.
 "S2":((-960,1800,2535),-1352),
 "S3":((-819,3185,5915),-3025),
 # First new repeated pair in the B=10000, 1001<=a stream.
 "U1":((-1071,1116,1134),-1054),
 "U2":((-1071,-1054,1116),1134),
 "V1":((-2178,2420,9075),-1470),
 "V2":((-1458,2268,7938),-2023),
 "W1":((-2254,-2162,2303),4900),
 "W2":((-2704,3042,8450),-2209),
 "X1":((-5929,-2541,9801),6069),
 "Y1":((-2835,5292,9450),-5290),
 "Z1":((-5408,-3840,7200),5415),
 "AA1":((-7488,-4680,7605),9800),
}

def projected_mask(p:int):
    base=set()
    path=ROOT/f"target_22224_direct_contact_deep13_boundary_p{p}.tsv"
    with path.open() as f:
        for r in csv.DictReader(f,delimiter="\t"):
            z=[int(r[k])%p for k in ("a","b","c","d")]
            base.add(tuple(sorted(x*x%p for x in z)))
    out=set()
    for key in base:
        for lam in range(1,p):
            out.add(tuple(sorted(lam*lam*x%p for x in key)))
    return base,out

rows=[]; lines=[]
for p in (11,13):
    base,mask=projected_mask(p)
    lines.append(f"p={p} base_square_multisets={len(base)} projected={len(mask)}")
    for name,(fixed,d0) in FIBERS.items():
        finite=[]
        for t in range(p):
            z=(*fixed,d0*t*t)
            key=tuple(sorted((x%p)**2%p for x in z))
            if key in mask: finite.append(t)
        infinity_key=tuple(sorted([0,0,0,(d0%p)**2%p]))
        infinity=(d0%p)!=0 and infinity_key in mask
        allzero=(d0%p)==0
        rows.append((name,p,",".join(map(str,finite)),int(infinity),int(allzero)))
        lines.append(f"{name} finite={finite} infinity={infinity} allzero_infinity_chart={allzero}")

out=ROOT/"target_22224_repeated_fibers_projected_contact.tsv"
with out.open("w") as f:
    f.write("fiber\tprime\tfinite_T_classes\tinfinity_allowed\tallzero_infinity_chart\n")
    for row in rows:f.write("\t".join(map(str,row))+"\n")
log=ROOT/"target_22224_repeated_fibers_projected_contact.log"
log.write_text("\n".join(lines)+"\n")
print(log.read_text(),end="")
