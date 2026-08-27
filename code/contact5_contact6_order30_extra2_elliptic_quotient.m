Q := Rationals();
P<u> := PolynomialRing(Q);
F := -3*(u-4)*(365*u^3 - 4044*u^2 + 14064*u - 14656);
C := HyperellipticCurve(F);
print "quotient quartic", F;
print "genus", Genus(C);
print "points_bound_100", Points(C : Bound := 100);
pt := C![4,0,1];
try
    E, mp := EllipticCurve(C, pt);
    print "E", E;
    print "map", mp;
    print "rank_bounds", RankBounds(E);
    T, phi := TorsionSubgroup(E);
    print "torsion", Invariants(T);
    print "generators", Generators(E);
    print "points_bound_100_on_E", Points(E : Bound := 100);
catch e
    print "elliptic conversion failed", e`Object;
end try;
quit;
