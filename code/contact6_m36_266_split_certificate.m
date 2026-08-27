//////////////////////////////////////////////////////////////////////
//  Exact split diagnostics for the known [2,6,6] examples.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

curves := [
    <"eps=1, r=21 and r=-1/3, a=187/21, b=-23/7",
     -111132*x^5 + 2646000*x^4 - 7101864*x^3 + 11226096*x^2 + 4630500*x>,
    <"eps=1, r=4/3, a=-1/42, b=-13/7",
     7112448*x^5 - 36091440*x^4 + 68732496*x^3 - 58231404*x^2 + 18522000*x>
];

print "Exact split diagnostics for known [2,6,6] examples";

Cprev := false;
for i in [1..#curves] do
    label := curves[i][1];
    f := curves[i][2];
    C := HyperellipticCurve(f);
    print "curve", i, label;
    print " factorization", Factorization(f);
    try
        print " automorphism_group", AutomorphismGroup(C);
    catch e
        print " automorphism_group unavailable", e`Object;
    end try;
    try
        subcovers := Degree2Subcovers(C);
        print " degree2_subcovers", #subcovers;
        for j in [1..#subcovers] do
            print "  subcover", j, subcovers[j][1], subcovers[j][2];
        end for;
    catch e
        print " degree2_subcovers unavailable", e`Object;
    end try;
    if i eq 2 then
        C1 := HyperellipticCurve(curves[1][2]);
        try
            print " is_isomorphic_to_curve1", IsIsomorphic(C1, C);
        catch e
            print " is_isomorphic_to_curve1 unavailable", e`Object;
        end try;
    end if;
end for;

quit;
