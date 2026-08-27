SetColumns(0); SetSeed(18);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
if not assigned log_file then
    log_file:="results/target_22224_B10000_a1000_newfiber_ranks.log";
end if;
SetLogFile(log_file:Overwrite:=true);

// Repeated square-ratio fibers that first appear after a=779 in the
// completed B=10000, 1<=a<=1000 all-sign stream.
fibers:=[
 <"S1",[Q!-4,9,30],Q!-2166/245,Q!7/19>,
 <"S2",[Q!-64,120,169],Q!-1352/15,Q!1/13>,
 <"S3",[Q!-9,35,65],Q!-3025/91,Q!5/11>
];

print "B10000_A1000_NEWFIBER_RANKS_START","fibers",#fibers;
for fam in fibers do
    name:=fam[1]; fixed:=fam[2]; d0:=fam[3]; t1:=fam[4];
    print "FIBER_START",name,"fixed",fixed,"d0",d0,"t1",t1;
    for ij in [<1,2>,<1,3>,<2,3>] do
        i,j:=Explode(ij);
        Ri:=(fixed[i]+d0*T^2)/(fixed[i]+d0);
        Rj:=(fixed[j]+d0*T^2)/(fixed[j]+d0);
        C:=HyperellipticCurve(Ri*Rj);
        P0:=C![Q!1,Q!1,Q!1]; EE,ph:=EllipticCurve(C,P0);
        E0,mm:=MinimalModel(EE); lo,hi:=RankBounds(E0);
        Tor:=TorsionSubgroup(E0);
        print "QUOTIENT",name,ij,"rank_bounds",lo,hi,
              "torsion",Invariants(Tor),"minimal",E0;
    end for;
    print "FIBER_DONE",name;
end for;
print "B10000_A1000_NEWFIBER_RANKS_DONE";
UnsetLogFile(); quit;
