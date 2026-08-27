Q := Rationals();
P<x> := PolynomialRing(Q);
f := -(3*x^2 - 12*x + 13)*
        (365*x^6 - 4044*x^5 + 18249*x^4 - 42664*x^3
         + 54039*x^2 - 34764*x + 8755);
C := HyperellipticCurve(f);
print "curve", C;
try
    print "Automorphisms", Automorphisms(C);
catch e
    print "Automorphisms failed", e`Object;
end try;
try
    print "Points bound 100", Points(C : Bound := 100);
catch e
    print "Points failed", e`Object;
end try;
try
    print "Chabauty no args", Chabauty(C);
catch e
    print "Chabauty failed", e`Object;
end try;
quit;
