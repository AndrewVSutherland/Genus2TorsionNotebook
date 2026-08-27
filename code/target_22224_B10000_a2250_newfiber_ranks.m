SetColumns(0); SetSeed(18);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
if not assigned log_file then
 log_file:="results/target_22224_B10000_a2250_newfiber_ranks.log";
end if;
SetLogFile(log_file:Overwrite:=true);
fibers:=[
 <"W1",[Q!-2254,-2162,2303],Q!4900,Q!7/5>,
 <"W2",[Q!-8,9,25],Q!-2209/338,Q!13/47>
];
print "B10000_A2250_NEWFIBER_RANKS_START","fibers",#fibers;
for fam in fibers do
 name:=fam[1]; fixed:=fam[2]; d0:=fam[3]; t1:=fam[4];
 print "FIBER_START",name,"fixed",fixed,"d0",d0,"t1",t1;
 for ij in [<1,2>,<1,3>,<2,3>] do
  i,j:=Explode(ij); Ri:=(fixed[i]+d0*T^2)/(fixed[i]+d0);
  Rj:=(fixed[j]+d0*T^2)/(fixed[j]+d0);
  C:=HyperellipticCurve(Ri*Rj); P0:=C![Q!1,Q!1,Q!1];
  EE,ph:=EllipticCurve(C,P0); E0,mm:=MinimalModel(EE);
  lo,hi:=RankBounds(E0); Tor:=TorsionSubgroup(E0);
  print "QUOTIENT",name,ij,"rank_bounds",lo,hi,
        "torsion",Invariants(Tor),"minimal",E0;
 end for;
 print "FIBER_DONE",name;
end for;
print "B10000_A2250_NEWFIBER_RANKS_DONE";
UnsetLogFile(); quit;
