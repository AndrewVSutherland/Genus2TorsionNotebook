//////////////////////////////////////////////////////////////////////
//  Contact-5/contact-5 symbolic obstruction.
//
//  We ask for a genus-2 odd quintic with two rational 5-contact
//  presentations at distinct rational x-coordinates.  After an affine
//  change in x, put the contacts at 0 and 1:
//
//      f = h0^2 - K*x^5 = h1^2 - K*(x-1)^5,
//
//  where h0,h1 are quadratics.  The leading coefficient of f forces
//  the same K in both presentations.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
R<A,B,C,e,d,c,K> := PolynomialRing(Q, 7, "grevlex");
PX<x> := PolynomialRing(R);

h0 := e*x^2 + d*x + c;
h1 := A*x^2 + B*x + C;
Phi := x^5 - (x-1)^5;
expr := h0^2 - h1^2 - K*Phi;
eqs := [Coefficient(expr, i) : i in [0..4]];

print "# normalized contact-5/contact-5 equations";
print "h0 =", h0;
print "h1 =", h1;
print "Phi = x^5-(x-1)^5 =", Phi;
for i in [1..#eqs] do
    printf "E%o = %o\n", i-1, eqs[i];
end for;

print "";
print "# sum/difference quotient";
u := h0 - h1;
v := h0 + h1;
print "u = h0-h1 =", u;
print "v = h0+h1 =", v;
print "u*v - K*Phi =", u*v - K*Phi;
print "Thus a nonzero solution gives a factorization K*Phi = u*v";
print "with deg(u),deg(v) <= 2.";

P<X> := PolynomialRing(Q);
PhiQ := X^5 - (X-1)^5;
print "";
print "# rational factorization check";
print "PhiQ =", PhiQ;
print "Factorization over Q:", Factorization(PhiQ);
print "IsIrreducible(PhiQ) =", IsIrreducible(PhiQ);
print "Discriminant(PhiQ) =", Discriminant(PhiQ);

K5<s5> := QuadraticField(5);
P5<Y> := PolynomialRing(K5);
Phi5 := Y^5 - (Y-1)^5;
print "";
print "# over Q(sqrt(5)) the obstruction splits into two quadratics";
print Factorization(Phi5);

Rr<r,Z> := PolynomialRing(Q, 2);
Pr<T> := PolynomialRing(Rr);
PhiR := T^5 - (T-r)^5;
print "";
print "# before normalizing the second contact to x=1";
print "x^5-(x-r)^5 =", PhiR;
print "factor =", Factorization(PhiR);
print "The r=0 component is the same-contact branch h1=+/-h0.";
print "For r != 0, substituting T=r*Z gives";
print Rr!(Evaluate(PhiR, r*Z)/r^5);
print "which is the same irreducible quartic in Z.";

assert IsIrreducible(PhiQ);
print "";
print "CONCLUSION: no nonzero rational distinct contact-5/contact-5";
print "presentation exists on an odd quintic.  K=0 is the square/boundary";
print "case, not a smooth genus-2 curve.";

quit;
