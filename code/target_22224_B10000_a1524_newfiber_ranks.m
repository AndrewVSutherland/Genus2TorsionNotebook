SetColumns(0); SetSeed(18);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
if not assigned log_file then
    log_file:="results/target_22224_B10000_a1524_newfiber_ranks.log";
end if;
SetLogFile(log_file:Overwrite:=true);
fibers:=[
 <"V1",[Q!-18,20,75],Q!-1470/121,Q!11/49>,
 <"V2",[Q!-9,14,49],Q!-2023/162,Q!3/17>
];
print "B10000_A1524_NEWFIBER_RANKS_START","fibers",#fibers;
for fam in fibers do
 name:=fam[1]; fixed:=fam[2]; d0:=fam[3]; t1:=fam[4];
 print "FIBER_START",name,"fixed",fixed,"d0",d0,"t1",t1;
 for ij in [<1,2>,<1,3>,<2,3>] do
  i,j:=Explode(ij);
  Ri:=(fixed[i]+d0*T^2)/(fixed[i]+d0);
  Rj:=(fixed[j]+d0*T^2)/(fixed[j]+d0);
  C:=HyperellipticCurve(Ri*Rj); P0:=C![Q!1,Q!1,Q!1];
  EE,ph:=EllipticCurve(C,P0); E0,mm:=MinimalModel(EE);
  lo,hi:=RankBounds(E0); Tor:=TorsionSubgroup(E0);
  print "QUOTIENT",name,ij,"rank_bounds",lo,hi,
        "torsion",Invariants(Tor),"minimal",E0;
 end for;
 print "FIBER_DONE",name;
end for;
print "B10000_A1524_NEWFIBER_RANKS_DONE";
UnsetLogFile(); quit;
