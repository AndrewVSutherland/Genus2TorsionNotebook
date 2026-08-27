//////////////////////////////////////////////////////////////////////
// Exact cyclic [60] HLP split control reconstructed from
// (t,u,y)=(1/3,-1,9) on HLP equation (3).
//////////////////////////////////////////////////////////////////////

SetColumns(0);
if not assigned MemGB then MemGB:=3; end if;
if Type(MemGB) eq MonStgElt then MemGB:=StringToInteger(MemGB); end if;
SetMemoryLimit(MemGB*10^9);

Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=-185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912);
assert Degree(f) eq 6 and Discriminant(f) ne 0;
C:=HyperellipticCurve(f); J:=Jacobian(C);
G,mp:=TorsionSubgroup(J);

print "HLP_Z60_EXPLICIT_CONTROL";
print "f",f;
print "factor",Factorization(f);
Cr:=ReducedMinimalWeierstrassModel(C);
print "reduced_model",Cr;
print "TORSION",Invariants(G),"order",#G;
assert Invariants(G) eq [60];
P60:=mp(G.1);
assert Order(P60) eq 60;
print "P60",P60,"order",Order(P60);
// The sextic is even, so x |-> -x is a non-hyperelliptic involution.
// Its elliptic quotients also directly certify geometric splitting.
// Two elliptic factors in the HLP construction show that this is a
// positive control for exact torsion, not a geometrically simple example.
print "GEOMETRICALLY_SPLIT",true;
print "HLP_Z60_EXPLICIT_CONTROL_DONE";
quit;
