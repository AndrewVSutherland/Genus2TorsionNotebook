SetColumns(0);
Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=183*(x^2+1)*(32*x^2+61*x+32)*(32*x^2-61*x+32);
print "CONTACT6_M612_HLP_SPLIT_SEED";
print "f",f;
print "factorization",Factorization(f);
C:=HyperellipticCurve(f); J:=Jacobian(C);
G,mp:=TorsionSubgroup(J);
print "torsion",Invariants(G),"order",#G;
for p in [5,7,11,13,17,19,23,29,31] do
 try
  fp:=ChangeRing(f,GF(p));
  if Discriminant(fp) eq 0 then continue; end if;
  Lp:=LPolynomial(ChangeRing(C,GF(p)));
  print "p",p,"L",Lp,"fac",Factorization(Lp);
 catch e
  continue;
 end try;
end for;
quit;
