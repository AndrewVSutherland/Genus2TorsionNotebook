//////////////////////////////////////////////////////////////////////
//  Halving the marked order-6 class on the contact-6 [3,6] family.
//
//  Put t=x-1 and
//      h6 = t^3 + A*t^2 + B*t + C,
//      f  = h6^2 - t^6.
//
//  The marked point is P=(t=0,y=C).  A rational half H of
//  D=P-infinity is represented by a quadratic Mumford divisor
//
//      u   = t^2 + p*t + q,
//      ell = r*t^2 + s*t - C,
//
//  satisfying
//
//      ell^2 - f = -2*A*t*u^2.
//
//  The sign ell(0)=-C makes the residual point (0,-C), hence
//  2H=P-infinity.  The opposite sign is equivalent after replacing H
//  by -H.
//
//  In the original contact-6 coordinates,
//      A=b+3, B=a+2*b+3, C=a+b+2,
//  so B=A+C-2.
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned output_file then
    output_file := "data/contact6_m36_halveD_symbolic.txt";
end if;

Q := Rationals();
R<A,C,p,q,r,s> := PolynomialRing(Q, 6);
P<t> := PolynomialRing(R);

B := A + C - 2;
h6 := t^3 + A*t^2 + B*t + C;
f := h6^2 - t^6;
u := t^2 + p*t + q;
ell := r*t^2 + s*t - C;

E := ell^2 - f + 2*A*t*u^2;
eqs := [Coefficient(E, i) : i in [1..4]];

// Two equations are linear in p and s on the open A*C != 0 chart.
p_expr := (A^2 + 2*B - r^2)/(4*A);
s_expr := A*q^2/C - B;

function PrimitiveNumerator(F)
    K := FieldOfFractions(R);
    return R!Numerator(K!F);
end function;

E2red := PrimitiveNumerator(Evaluate(eqs[2], [A,C,p_expr,q,r,s_expr]));
E3red := PrimitiveNumerator(Evaluate(eqs[3], [A,C,p_expr,q,r,s_expr]));

print "Contact-6 marked-class halving equations";
print "h6 =", h6;
print "f =", f;
print "u =", u;
print "ell =", ell;
print "";
print "Equations H1..H4:";
for i in [1..#eqs] do
    print Sprintf("H%o =", i), eqs[i];
end for;
print "";
print "Open-chart substitutions:";
print "p =", p_expr;
print "s =", s_expr;
print "";
print "Remaining equations after eliminating p,s:";
print "K2 =", E2red;
print "K3 =", E3red;
print "";
print "Nonboundary conditions:";
print "A*C*q != 0, discriminant(f) != 0, gcd(u,f)=1";

out := Open(output_file, "w");
fprintf out, "# Contact-6 marked-class halving equations\n";
fprintf out, "h6 = %o\n", h6;
fprintf out, "f = %o\n", f;
fprintf out, "u = %o\n", u;
fprintf out, "ell = %o\n\n", ell;
for i in [1..#eqs] do
    fprintf out, "H%o = %o\n", i, eqs[i];
end for;
fprintf out, "\nOpen chart A*C != 0:\n";
fprintf out, "p = %o\n", p_expr;
fprintf out, "s = %o\n", s_expr;
fprintf out, "K2 = %o\n", E2red;
fprintf out, "K3 = %o\n", E3red;
fprintf out, "\nA=b+3, C=a+b+2, B=A+C-2\n";
fprintf out, "Nonboundary: A*C*q != 0, discriminant(f) != 0, gcd(u,f)=1\n";
delete out;

print "Wrote", output_file;
quit;
