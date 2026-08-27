//////////////////////////////////////////////////////////////////////
//  Exact split certificate for the M(2,12)+extra-3 hit.
//
//  The hit is the curve
//
//      y^2 = 5668704*x^5 - 22143375*x^4 + 36098622*x^3
//            - 30305259*x^2 + 12990780*x - 2259900.
//
//  Magma's Degree2Subcovers produces two degree-2 maps to elliptic curves.
//  The raw elliptic curves minimize to Cremona labels 90c3 and 510g1.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
P<X> := PolynomialRing(Q);

f := 5668704*X^5 - 22143375*X^4 + 36098622*X^3
     - 30305259*X^2 + 12990780*X - 2259900;
C := HyperellipticCurve(f);

print "M(2,12)+extra3 split certificate";
print "C:", C;
print "AutomorphismGroup:", AutomorphismGroup(C);

subcovers := Degree2Subcovers(C);
print "degree2_subcovers", #subcovers;
for i in [1..#subcovers] do
    E := subcovers[i][1];
    phi := subcovers[i][2];
    print "subcover", i;
    print "  E_raw =", E;
    print "  map =", phi;
end for;

print "Minimal models identified in Sage:";
print "  subcover 1 -> Cremona 90c3:";
print "    y^2 + x*y + y = x^3 - x^2 - 122*x + 1721";
print "  subcover 2 -> Cremona 510g1:";
print "    y^2 + x*y = x^3 + 25*x - 375";

quit;
