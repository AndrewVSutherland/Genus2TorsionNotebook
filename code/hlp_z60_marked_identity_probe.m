//////////////////////////////////////////////////////////////////////
// Probe the exact order-3, order-4, and order-5 multiples of the
// HLP cyclic-[60] control generator.  This is intentionally compact and
// deterministic; the final self-contained verifier is written separately.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
SetMemoryLimit(2*10^9);
Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=-185*(125*x^2-1728)*(2000*x^4-48525*x^2+294912);
C:=HyperellipticCurve(f); J:=Jacobian(C);
G,mp:=TorsionSubgroup(J);
assert Invariants(G) eq [60];
P60:=J![x^2+74/15*x+128/25,-1081510/9*x-1226624/3]; assert Order(P60) eq 60;
for pair in [<"D5",12>,<"D3",20>,<"D4",15>] do
    D:=pair[2]*P60;
    print pair[1],"order",Order(D),"D",D;
    print "Eltseq",Eltseq(D);
end for;
print "D2",30*P60,"Eltseq",Eltseq(30*P60);
quit;
