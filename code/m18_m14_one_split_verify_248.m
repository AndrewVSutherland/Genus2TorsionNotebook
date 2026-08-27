SetColumns(0);

Q := Rationals();
P<x> := PolynomialRing(Q);

f := 7061463847622250*x^5
     + 104632219276049025*x^4
     + 135735215960638800*x^3
     + 188573481843278400*x^2
     + 51200550567936000*x;

C := HyperellipticCurve(f);
J := Jacobian(C);
G, phi := TorsionSubgroup(J);
print "torsion", Invariants(G);

for p in [47,3,5,7,11,13,17,19,23,29,31,37,41,43,53,59,61,67,71,73] do
    try
        fp := ChangeRing(f, GF(p));
        if Degree(fp) eq 5 and Discriminant(fp) ne 0 then
            Lp := LPolynomial(ChangeRing(C, GF(p)));
            print "p", p, "L", Lp, "fac", Factorization(Lp);
        end if;
    catch e
        dummy := 0;
    end try;
end for;

try
    Cr := ReducedMinimalWeierstrassModel(C);
    print "ReducedMinimalWeierstrassModel", Cr;
catch e
    print "ReducedMinimalWeierstrassModel error", e`Object;
end try;

try
    Cr2 := MinimalWeierstrassModel(C);
    print "MinimalWeierstrassModel", Cr2;
catch e
    print "MinimalWeierstrassModel error", e`Object;
end try;

try
    Cr3 := SimplifiedModel(C);
    print "SimplifiedModel", Cr3;
catch e
    print "SimplifiedModel error", e`Object;
end try;

quit;
