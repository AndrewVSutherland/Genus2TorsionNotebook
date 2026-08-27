// twists.m — Lane D task 3: exact torsion of quadratic twists of the
// (2,2,2,12) hit curve, squarefree D with |D| <= 30.  Full 2-torsion is
// twist-invariant here (all 6 Weierstrass points rational), so every twist
// has torsion >= (2,2,2,2); the question is whether any keeps the 4 or 3 part.
P<x> := PolynomialRing(Rationals());
A := [1,1,1,2,2];
B := [282322361376, -8243383980, -64241207724, -114724491840, 561915878400];
f0 := &*[A[i] + B[i]*x : i in [1..5]];
C0 := HyperellipticCurve(f0);
for D in [d : d in [-30..30] | d ne 0 and d ne 1 and IsSquarefree(d)] do
    CD := QuadraticTwist(C0, D);
    CR := CD;
    try
        CR := ReducedMinimalWeierstrassModel(CD);
    catch e
        CR := CD;
    end try;
    CR := IntegralModel(SimplifiedModel(CR));
    inv := Invariants(TorsionSubgroup(Jacobian(CR)));
    printf "D=%o: torsion %o order %o\n", D, inv, &*inv;
end for;
print "TWISTS DONE";
quit;
