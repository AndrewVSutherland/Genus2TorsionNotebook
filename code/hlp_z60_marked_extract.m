//////////////////////////////////////////////////////////////////////
// Extract the order 5, 3, and 4 layers from the exact HLP Z/60 seed.
// This is a small diagnostic used by the transverse-deformation analysis.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);
Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=-185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912);
C:=HyperellipticCurve(f); J:=Jacobian(C);
G,mp:=TorsionSubgroup(J);
assert Invariants(G) eq [60];
T:=mp(G.1);

print "HLP_Z60_MARKED_EXTRACT";
print "F",f;
print "T60",T,"order",Order(T);
for n in [12,20,15,30] do
    D:=n*T;
    print "MULTIPLE",n,"class",D,"order",Order(D);
    print "ELTSEQ",Eltseq(D);
end for;
print "DONE";
quit;
