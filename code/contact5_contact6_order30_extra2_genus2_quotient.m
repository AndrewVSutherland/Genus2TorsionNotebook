//////////////////////////////////////////////////////////////////////
//  Genus-2 quotient of the extra-2 genus-3 curve by sigma*hyperelliptic.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();
P<u> := PolynomialRing(Q);

H4 := (u-4)*(365*u^3 - 4044*u^2 + 14064*u - 14656);
D := 9*u^2 - 60*u + 84;
F := -3*D*H4;
C := HyperellipticCurve(F);
J := Jacobian(C);

print "genus2 quotient", F;
print "factorization", Factorization(F);
print "genus", Genus(C);
print "points bound 1000", Points(C : Bound := 1000);
try
    print "rank_bounds", RankBounds(J);
catch e
    print "rank_bounds failed", e`Object;
end try;
try
    T, phi := TorsionSubgroup(J);
    print "torsion", Invariants(T);
catch e
    print "torsion failed", e`Object;
end try;

// Images of the known genus-3 point R=1,Y=16.
// u=2; quotient Z^2 = F(2).
print "F(2)", Evaluate(F, 2);

try
    D0 := J!0;
    pts_chab := Chabauty(D0 : ptC := C![2,0,1]);
    print "chabauty_zero", pts_chab;
catch e
    print "chabauty_zero failed", e`Object;
end try;
try
    pts_chab2 := Chabauty(C![2,0,1]);
    print "chabauty_point", pts_chab2;
catch e
    print "chabauty_point failed", e`Object;
end try;

quit;
