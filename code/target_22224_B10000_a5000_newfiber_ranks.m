SetColumns(0); SetSeed(18);
Q:=Rationals(); R<T>:=PolynomialRing(Q);
if not assigned log_file then log_file:="results/target_22224_B10000_a5000_newfiber_ranks.log"; end if;
SetLogFile(log_file:Overwrite:=true);
fixed:=[Q!-64,-40,65]; d0:=Q!9800/117; t1:=Q!9/7;
print "FIBER_START","AA1","fixed",fixed,"d0",d0,"t1",t1;
for ij in [<1,2>,<1,3>,<2,3>] do
 i,j:=Explode(ij); Ri:=(fixed[i]+d0*T^2)/(fixed[i]+d0);
 Rj:=(fixed[j]+d0*T^2)/(fixed[j]+d0);
 C:=HyperellipticCurve(Ri*Rj); P0:=C![Q!1,Q!1,Q!1];
 EE,ph:=EllipticCurve(C,P0); E0,mm:=MinimalModel(EE);
 lo,hi:=RankBounds(E0); Tor:=TorsionSubgroup(E0);
 print "QUOTIENT","AA1",ij,"rank_bounds",lo,hi,
       "torsion",Invariants(Tor),"minimal",E0;
end for;
print "FIBER_DONE","AA1"; UnsetLogFile(); quit;
