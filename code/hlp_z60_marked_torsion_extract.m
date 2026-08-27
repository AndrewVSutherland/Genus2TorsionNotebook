//////////////////////////////////////////////////////////////////////
// Extract the cyclic order-60 generator and its marked primary parts
// from the exact split HLP control.  These data seed the simultaneous
// order-5/order-3/order-4 deformation calculation.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=2; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=-185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912);
C:=HyperellipticCurve(f); J:=Jacobian(C);
G,mp:=TorsionSubgroup(J);
assert Invariants(G) eq [60];
P60:=J![x^2+74/15*x+128/25,-1081510/9*x-1226624/3];
assert Order(P60) eq 60;

print "HLP_Z60_MARKED_TORSION_EXTRACT";
print "f",f;
print "P60",P60,"order",Order(P60);
for n in [12,20,15,30] do
    D:=n*P60;
    print "MULT",n,"order",Order(D),"jacobian",D;
end for;
print "HLP_Z60_MARKED_TORSION_EXTRACT_DONE";
quit;
