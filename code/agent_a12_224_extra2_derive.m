
//////////////////////////////////////////////////////////////////////
//  Z/2 x Z/24 extra-2 locus on the A(12) chart.
//
//  f = R*F (deg 6), R quadratic, F = Q^2 + R*ell^2 quartic; generic
//  factor type [2,4] => 2-rank 1.  For the target [2,24] we need
//  2-rank 2, i.e. factor type [2,2,2]: R irreducible, F split into two
//  rational quadratics (F -> [2,2]).  This adds one independent rational
//  2-torsion class T_F, so the halving funnel can target [2,24] rather
//  than [24].
//
//  This script:
//   1. builds R, F over Q(p,z,r), confirms deg/generic factor type;
//   2. forms the "F splits into two rational quadratics" condition via
//      its cubic resolvent (a root that is a perfect square), and
//      factors it over Q(p,z,r) to find a clean parametrization;
//   3. also reports disc(F) and the R-split factor for reference.
//////////////////////////////////////////////////////////////////////

SetColumns(0);
Q := Rationals();

K3<p,z,r> := RationalFunctionField(Q, 3);
PX<x> := PolynomialRing(K3);

s := (z^2 - 4*p^2 + 1)/(2*z);
t := (z^2 + 4*p^2 - 1)^2/(8*p^2*z);
mu := ((s^2 - 1)*(2*p*r + 1) - p^2*(2*s*t - 4))/(4*p^3);
lambda := (4 - mu^2)*p^2/(s^2 - 1);
T1 := p*x + r;
Rpol := (T1^2 + x - 1)/lambda;
ell := s*x + t;
Qpol := 2*T1 + mu*Rpol;
F := Rpol*x^2 + 4*(Rpol + x - 1)*(Rpol - 1);
assert F eq Qpol^2 + Rpol*ell^2;

printf "deg R = %o, deg F = %o\n", Degree(Rpol), Degree(F);

// normalize F to monic quartic x^4 + c3 x^3 + c2 x^2 + c1 x + c0
lcF := LeadingCoefficient(F);
Fm := F/lcF;
c3 := Coefficient(Fm,3); c2 := Coefficient(Fm,2);
c1 := Coefficient(Fm,1); c0 := Coefficient(Fm,0);
printf "F monic coeffs c3,c2,c1,c0 computed (rational functions of p,z,r)\n";

// depress: x = y - c3/4; y^4 + P y^2 + Qq y + Rr
Pp := c2 - 3*c3^2/8;
Qq := c1 - c2*c3/2 + c3^3/8;
Rr := c0 - c1*c3/4 + c2*c3^2/16 - 3*c3^4/256;

// F splits into two rational quadratics <=> the resolvent cubic
//   y^3 + 2P y^2 + (P^2 - 4R) y - Q^2 = 0   (y = beta^2 form)
// has a rational root that is a perfect square, OR (cleaner) parametrize
// directly: F = (x^2 + a1 x + b1)(x^2 + a2 x + b2) with a1 + a2 = c3,
// a1 a2 + b1 + b2 = c2, a1 b2 + a2 b1 = c1, b1 b2 = c0.
// Use the standard resolvent whose rational root gives the split:
resolv := Pp^2 - 4*Rr;   // as coefficients; build cubic in u
Pu<u> := PolynomialRing(K3);
res := u^3 + 2*Pp*u^2 + (Pp^2 - 4*Rr)*u - Qq^2;
print "resolvent cubic (u = square-of-half-difference param):";
printf "  res = u^3 + (%o) u^2 + ... (deg-3 in u over Q(p,z,r))\n", 2*Pp;

// The split condition is: res has a rational root.  Compute its
// discriminant / try to factor res over Q(p,z,r)[u].
fac := Factorization(res);
printf "resolvent factors over Q(p,z,r): %o\n",
    [<Degree(ff[1]), ff[2]> : ff in fac];
for ff in fac do
    if Degree(ff[1]) eq 1 then
        // rational root u0 = -constant/leading
        u0 := -Coefficient(ff[1],0)/Coefficient(ff[1],1);
        printf "  RATIONAL RESOLVENT ROOT u0 (F may split for free!):\n";
        num := Numerator(u0); den := Denominator(u0);
        printf "   u0 num factors: %o\n", Factorization(num);
        printf "   u0 den factors: %o\n", Factorization(den);
    end if;
end for;

// direct approach: try to factor F itself over Q(p,z,r)
print "\ndirect factorization of F over Q(p,z,r):";
facF := Factorization(F);
printf "F factors: %o\n", [<Degree(ff[1]), ff[2]> : ff in facF];

// disc(F) for reference
dF := Discriminant(Fm);
printf "\ndisc(F) squarefree-relevant factors:\n";
ndF := Numerator(dF); ddF := Denominator(dF);
print "  num:", Factorization(ndF);
print "  den:", Factorization(ddF);

print "DONE";
quit;
