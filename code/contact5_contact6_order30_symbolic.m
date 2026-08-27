//////////////////////////////////////////////////////////////////////
//  Simultaneous contact-5/contact-6 construction for order 30.
//
//  Seek an odd quintic f with
//
//      f = h6^2 - (x-1)^6 = h5^2 - K*x^5,
//
//  where
//
//      h6 = x^3 + A*x^2 + B*x + C,
//      h5 = e*x^2 + d*x + c.
//
//  Then y-h6 has divisor 6*(1,h6(1))-6*infinity, and y-h5 has
//  divisor 5*(0,c)-5*infinity.  Smooth nonboundary specializations
//  should therefore have a rational point of order lcm(5,6)=30,
//  unless the two classes are accidentally dependent.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/contact5_contact6_order30_symbolic.txt";
end if;
if not assigned do_groebner then
    do_groebner := false;
end if;

Q := Rationals();
R<A,B,C,e,d,c,K> := PolynomialRing(Q, 7, "grevlex");
PX<x> := PolynomialRing(R);

h6 := x^3 + A*x^2 + B*x + C;
h5 := e*x^2 + d*x + c;
f6 := h6^2 - (x-1)^6;
expr := f6 - h5^2 + K*x^5;
eqs := [Coefficient(expr, i) : i in [0..5]];

out := Open(output_file, "w");
fprintf out, "# Simultaneous contact-5/contact-6 equations\n\n";
fprintf out, "h6 = %o\n", h6;
fprintf out, "h5 = %o\n", h5;
fprintf out, "f = h6^2 - (x-1)^6\n\n";
for i in [1..#eqs] do
    fprintf out, "E%o = %o\n", i-1, eqs[i];
end for;
fprintf out, "\nTop equation solves K = %o\n", -Coefficient(f6,5);

// Parametrize C^2-c^2=1 by u:
//      C=(u^2+1)/(2u), c=(u^2-1)/(2u).
// On c != 0, E1 solves d=(B*C+3)/c.
// On E4, solve B=(e^2-A^2+15)/2.
Ru<u,A1,e1> := PolynomialRing(Q, 3, "grevlex");
Ku := FieldOfFractions(Ru);
PXu<X> := PolynomialRing(Ku);

Cu := (u^2 + 1)/(2*u);
cu := (u^2 - 1)/(2*u);
Bu := (e1^2 - A1^2 + 15)/2;
du := (Bu*Cu + 3)/cu;
h6u := X^3 + A1*X^2 + Bu*X + Cu;
h5u := e1*X^2 + du*X + cu;
fu := h6u^2 - (X-1)^6;
Gu := fu - h5u^2 - Coefficient(fu,5)*X^5;
F2 := Numerator(Coefficient(Gu, 2));
F3 := Numerator(Coefficient(Gu, 3));
F2 := Ru!F2;
F3 := Ru!F3;

fprintf out, "\nParametrized open chart C^2-c^2=1, c != 0:\n";
fprintf out, "C = (u^2+1)/(2u)\n";
fprintf out, "c = (u^2-1)/(2u)\n";
fprintf out, "B = (e^2-A^2+15)/2\n";
fprintf out, "d = (B*C+3)/c\n\n";
fprintf out, "Remaining equations in (u,A,e):\n";
fprintf out, "F2 = %o\n", F2;
fprintf out, "F3 = %o\n", F3;
fprintf out, "F2 factorization:\n";
for pair in Factorization(F2) do
    fprintf out, "  exponent %o: %o\n", pair[2], pair[1];
end for;
fprintf out, "F3 factorization:\n";
for pair in Factorization(F3) do
    fprintf out, "  exponent %o: %o\n", pair[2], pair[1];
end for;
fprintf out, "gcd(F2,F3) = %o\n", GCD(F2,F3);

res_e := Resultant(F2, F3, e1);
fprintf out, "\nResultant eliminating e: degree %o, total degree %o, terms %o\n",
        Degree(res_e), TotalDegree(res_e), #Terms(res_e);
fprintf out, "factorization of Res_e(F2,F3):\n";
for pair in Factorization(res_e) do
    fprintf out, "  exponent %o degree %o total_degree %o terms %o: %o\n",
            pair[2], Degree(pair[1]), TotalDegree(pair[1]), #Terms(pair[1]), pair[1];
end for;

if do_groebner then
    I := ideal<Ru | F2, F3>;
    fprintf out, "\nComputing Groebner basis of <F2,F3> in Ru...\n";
    time G := GroebnerBasis(I);
    fprintf out, "GB length %o\n", #G;
    for i in [1..#G] do
        fprintf out, "GB%o = %o\n", i, G[i];
    end for;
end if;

delete out;
print "Wrote", output_file;
quit;
