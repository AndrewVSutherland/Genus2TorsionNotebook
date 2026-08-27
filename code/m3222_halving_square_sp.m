//////////////////////////////////////////////////////////////////////
//  Symmetric-parameter square-quartic equations for halving Q in the
//  odd M_1(8,2,2) family.
//
//  Write s=u+v and p=uv.  If H=(a,b) is a degree-2 reduced divisor
//  whose double is Q=(-1, p*(s+1)), then
//
//      f_{s,p}(X) - L*(X+1)*(X^2 + A*X + B)^2
//
//  must be a square of a quadratic.  The leading coefficient gives
//      L = coeff_X^5(f_{s,p}) = -(p-s+1)*(s+2).
//
//  This script forms the residual quartic coefficients r0..r4 and the
//  universal square-quartic equations in s,p,A,B.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
R<s,p,A,B> := PolynomialRing(Q, 4, "grevlex");
PX<X> := PolynomialRing(R);

qtilde := -X^2 + (p*s - s^2 + 2*p - s - 2)*X - (s^2 - p + s + 1);
f := ((p-s+1)*X^2 + (2-s)*X + 1)*((s+2)*X + 1)*qtilde;
a := X^2 + A*X + B;
L := Coefficient(f, 5);
res := f - L*(X+1)*a^2;
r := [Coefficient(res, i) : i in [0..4]];

F := [
    8*r[1]^2*r[4] - 4*r[1]*r[2]*r[3] + r[2]^3,
    16*r[1]^2*r[5] + 2*r[1]*r[2]*r[4] - 4*r[1]*r[3]^2 + r[2]^2*r[3],
    8*r[1]*r[2]*r[5] - 4*r[1]*r[3]*r[4] + r[2]^2*r[4],
    r[1]*r[4]^2 - r[2]^2*r[5],
    8*r[1]*r[4]*r[5] - 4*r[2]*r[3]*r[5] + r[2]*r[4]^2,
    16*r[1]*r[5]^2 + 2*r[2]*r[4]*r[5] - 4*r[3]^2*r[5] + r[3]*r[4]^2,
    8*r[2]*r[5]^2 - 4*r[3]*r[4]*r[5] + r[4]^3
];

print "Symmetric halving square-quartic setup";
print "variables: s,p,A,B";
print "f", f;
print "L", L;
print "residual coefficients r0..r4";
for i in [1..5] do
    print "r", i-1, r[i];
end for;

print "square equations degrees", [TotalDegree(g) : g in F];
for i in [1..#F] do
    fac := Factorization(F[i]);
    print "F", i, "degree", TotalDegree(F[i]), "terms", #Terms(F[i]), "factor_count", #fac;
    print fac;
end for;

if assigned do_groebner and do_groebner then
    I := ideal<R | F>;
    print "Computing Groebner basis of square-quartic ideal...";
    time G := GroebnerBasis(I);
    print "GB length", #G;
    for i in [1..#G] do
        print "GB", i, "degree", TotalDegree(G[i]), "terms", #Terms(G[i]);
        print G[i];
    end for;
end if;

quit;
