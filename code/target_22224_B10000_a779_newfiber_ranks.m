SetColumns(0); SetSeed(18);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
if not assigned log_file then log_file:="results/target_22224_B10000_a779_newfiber_ranks.log"; end if;
SetLogFile(log_file:Overwrite:=true);
fibers:=[
 <"A",[Q!-1,2,25],Q!-529/338,Q!13/23>,
 <"B",[Q!-44,-16,49],Q!676/11,Q!43/13>,
 <"C",[Q!-17,-16,50],Q!338/17,Q!53/13>,
 <"D",[Q!-7,16,56],Q!-392/25,Q!5/13>,
 <"E",[Q!-23,-18,64],Q!578/23,Q!65/17>,
 <"F",[Q!-3,5,75],Q!-405/121,Q!11/21>,
 <"G",[Q!-47,-4,98],Q!2738/47,Q!65/37>
];
print "B10000_A779_NEWFIBER_RANKS_START","fibers",#fibers;
for fam in fibers do name:=fam[1];fixed:=fam[2];d0:=fam[3];t1:=fam[4];
 print "FIBER_START",name,"fixed",fixed,"d0",d0,"t1",t1;
 for ij in [<1,2>,<1,3>,<2,3>] do i,j:=Explode(ij);
  Ri:=(fixed[i]+d0*T^2)/(fixed[i]+d0);Rj:=(fixed[j]+d0*T^2)/(fixed[j]+d0);
  C:=HyperellipticCurve(Ri*Rj);P0:=C![Q!1,Q!1,Q!1];EE,ph:=EllipticCurve(C,P0);
  E0,mm:=MinimalModel(EE);lo,hi:=RankBounds(E0);Tor:=TorsionSubgroup(E0);
  print "QUOTIENT",name,ij,"rank_bounds",lo,hi,"torsion",Invariants(Tor),"minimal",E0;
 end for; print "FIBER_DONE",name;
end for;
print "B10000_A779_NEWFIBER_RANKS_DONE";UnsetLogFile();quit;
