//////////////////////////////////////////////////////////////////////
// Extract Mumford representatives for the exact [6,12] torsion on the
// explicit HLP split control.  The order-3 and order-4 multiples provide
// marked data for a possible deformation away from the split locus.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q:=Rationals(); P<x>:=PolynomialRing(Q);
f:=183*(x^2+1)*(32*x^2+61*x+32)*(32*x^2-61*x+32);
C:=HyperellipticCurve(f); J:=Jacobian(C);
G,mp:=TorsionSubgroup(J);
print "M612_HLP_MARKED_TORSION_EXTRACT";
print "invariants",Invariants(G),"f",f;
gens:=OrderedGenerators(G);
for i in [1..#gens] do
    g:=gens[i]; D:=mp(g);
    print "GEN",i,"abstract",g,"order",Order(g),"jacobian",D;
    for n in [2,3,4,6] do
        print " MULT",n,"order",Order(n*g),"jacobian",mp(n*g);
    end for;
end for;

// Enumerate the four nonzero order-3 elements up to sign and every
// order-4 element, recording their reduced Mumford data as printed by Magma.
for g in G do
    if Order(g) eq 3 then
        print "ORDER3",g,"jacobian",mp(g);
    end if;
end for;
for g in G do
    if Order(g) eq 4 then print "ORDER4",g,"jacobian",mp(g); end if;
end for;
quit;
