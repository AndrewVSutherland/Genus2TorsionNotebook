//////////////////////////////////////////////////////////////////////
//  Fixed-slice elimination for the contact-6 extra-root route.
//
//  Default slice: eps=1, r=2.  Set fixed_eps/fixed_r_num/fixed_r_den
//  on the command line to analyze another slice.
//
//      magma -b code/contact6_m36_extra_root_slice_component.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

if not assigned fixed_eps then
    fixed_eps := 1;
elif Type(fixed_eps) eq MonStgElt then
    fixed_eps := StringToInteger(fixed_eps);
end if;
if not assigned fixed_r_num then
    fixed_r_num := 2;
elif Type(fixed_r_num) eq MonStgElt then
    fixed_r_num := StringToInteger(fixed_r_num);
end if;
if not assigned fixed_r_den then
    fixed_r_den := 1;
elif Type(fixed_r_den) eq MonStgElt then
    fixed_r_den := StringToInteger(fixed_r_den);
end if;
if not assigned output_file then
    output_file := "data/contact6_m36_extra_root_slice_component.txt";
end if;

Q := Rationals();
Z := Integers();
R<b,M,U,v> := PolynomialRing(Q, 4);

r := Q!fixed_r_num / Q!fixed_r_den;
eps := Q!fixed_eps;
a := (eps*(r-1)^3 - 1 - r^3 - b*r^2)/r;

c1 := 2*a + 6;
c2 := a^2 + 2*b - 15;
c3 := 2*a*b + 22;
c4 := 2*a + b^2 - 15;
c5 := 2*b + 6;

B := c5*M + 3*U;
Delta := 4*c4*M + 12*(U^2 + v^2) - B^2;
F3 := B*Delta + 16*v^3 - 8*c3*M - 8*U^3 - 48*U*v^2;
F2 := Delta^2 + 64*B*v^3 - 64*c2*M - 192*(U^2*v^2 + v^4);
F1 := Delta*v^3 - 4*c1*M - 12*U*v^4;

function PrimitivePolynomial(f)
    den := LCM([Denominator(c) : c in Coefficients(f)]);
    g := Parent(f)!(den*f);
    nums := [Z!c : c in Coefficients(g)];
    cont := GCD([Abs(n) : n in nums | n ne 0]);
    if cont gt 1 then
        g := Parent(f)!(g/cont);
    end if;
    return g;
end function;

Fs := [PrimitivePolynomial(R!F) : F in [F1,F2,F3]];
S<M0,U0,v0> := PolynomialRing(Q, 3);
phi := hom<R -> S | 0, M0, U0, v0>;
Res12 := PrimitivePolynomial(phi(Resultant(Fs[1], Fs[2], b)));
Res13 := PrimitivePolynomial(phi(Resultant(Fs[1], Fs[3], b)));
common := GCD(Res12, Res13);

out := Open(output_file, "w");
fprintf out, "Extra-root fixed-slice component check\n";
fprintf out, "eps = %o\nr = %o\n\n", eps, r;
fprintf out, "F degrees = %o\n", [TotalDegree(f) : f in Fs];
fprintf out, "F terms = %o\n\n", [#Terms(f) : f in Fs];
fprintf out, "Res(F1,F2;b): degree %o, terms %o, factor summary %o\n",
        TotalDegree(Res12), #Terms(Res12),
        [<TotalDegree(fe[1]), #Terms(fe[1]), fe[2]> : fe in Factorization(Res12)];
fprintf out, "Res(F1,F3;b): degree %o, terms %o, factor summary %o\n",
        TotalDegree(Res13), #Terms(Res13),
        [<TotalDegree(fe[1]), #Terms(fe[1]), fe[2]> : fe in Factorization(Res13)];
fprintf out, "gcd degree %o, terms %o, factor summary %o\n",
        TotalDegree(common), #Terms(common),
        [<TotalDegree(fe[1]), #Terms(fe[1]), fe[2]> : fe in Factorization(common)];
fprintf out, "gcd = %o\n", common;
delete out;

print "Extra-root fixed-slice component check";
print "eps", eps, "r", r;
print "F degrees", [TotalDegree(f) : f in Fs], "terms", [#Terms(f) : f in Fs];
print "Res12 degree", TotalDegree(Res12), "terms", #Terms(Res12),
      "factor_summary", [<TotalDegree(fe[1]), #Terms(fe[1]), fe[2]> : fe in Factorization(Res12)];
print "Res13 degree", TotalDegree(Res13), "terms", #Terms(Res13),
      "factor_summary", [<TotalDegree(fe[1]), #Terms(fe[1]), fe[2]> : fe in Factorization(Res13)];
print "gcd degree", TotalDegree(common), "terms", #Terms(common),
      "factor_summary", [<TotalDegree(fe[1]), #Terms(fe[1]), fe[2]> : fe in Factorization(common)];
print "gcd", common;
print "wrote", output_file;

quit;
