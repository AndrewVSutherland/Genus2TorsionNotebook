// claude_ov_88chi_symbolic.m -- Lane 5 (overnight 2026-07-25).
// THE SYMBOLIC LIFT LAYER for [8,8] on Nicholls' Lambda_334 family, over Q(s,t,v).
//
// Verifies, identically in (s,t,v):
//   (1) b is a perfect square:  A*N = (s^2(Av+1)^2 + t^2(1-t^2))^2, b = b0^2 with
//       b0 = A t (s^2 A v^2 - 1)/Delta,  Delta := s^2(Av+1)^2 + t^2(1-t^2).
//       Hence l2 = -x^2+bc splits and the sextic has factor type [1,1,2,2].
//   (2) Res(L1,L2) = (bc-a)^2 (b-c)^2 / e2^2      (a square)
//       Res(L2,L3) = (bc-a)^2                     (a square)
//       L1(r)  = (bc-a)(b0-t)^2/e2,  L1(-r) = (bc-a)(b0+t)^2/e2,  r := b0*t = sqrt(bc)
//   (3) c - b  ==  Psi := v(1+s^2 v)(1+A v)   modulo squares.
//       Psi = square is EXACTLY the double-stage-1 condition (the class {c,oo} being
//       2-divisible on J2), i.e. v = 1/z for a rational point (z,w) of the elliptic
//       fibration E_{s,t} : w^2 = z(z+s^2)(z+A).
//   (4) the order-4 Mumford u-polynomials above T2 = [L2] have the shape
//         u(x) = x^2 + w x - (w t^2 + b c)
//       and the halving system  (A+Bx)^2 L2 - d2 L1 L3 = k u^2  forces
//         k = -d2 (bc-a)^2 / (e2 t^2 w^2).
//
// Combined with 2-rank(J1) = 2 identically, this proves
//   [8,8] on Lambda_334  ==>  Psi is a square,
// i.e. the [8,8] locus is contained in the double-stage-1 elliptic fibration.
SetColumns(0);
SetMemoryLimit(20*10^9);
Qq := Rationals();
R<s,t,v> := PolynomialRing(Qq, 3);
F := FieldOfFractions(R);
Px<x> := PolynomialRing(F);

A  := s^2 - t^4 + t^2;
uu := (-s^2*A*v^2 - 2*A*v - 1)/(-s^2*t*A*v^2 + t);
a  := A/(1 - t^2);
b  := A/(uu^2*s^2 + 1 - t^2);
c  := t^2;
d2 := A*(s^2*uu^2 + t^4 - 2*t^2 + 1)
        *(s^4*uu^2 - s^2*t^2*uu^2 + s^2*uu^2 - t^6 + 3*t^4 - 3*t^2 + 1);

e2 := -a + b + c - 1;
e1 := 2*a - 2*b*c;
e0 := a*b*c - a*b - a*c + b*c;
l1 := e2*x^2 + e1*x + e0;
l2 := -x^2 + b*c;
l3 := x^2 - a;
L1 := x^2 + (e1/e2)*x + (e0/e2);
L2 := x^2 - b*c;
L3 := x^2 - a;

printf "=== (1) b is a perfect square ===\n";
Delta := s^2*(A*v+1)^2 + t^2*(1-t^2);
b0 := A*t*(s^2*A*v^2 - 1)/Delta;
printf "b - b0^2 = %o\n", b - b0^2;
r := b0*t;                      // r^2 = b*c
printf "r^2 - b*c = %o\n", r^2 - b*c;

printf "=== (2) resultants and L1(+-r) ===\n";
printf "Res(L1,L2) - (b*c-a)^2*(b-c)^2/e2^2 = %o\n",
       Resultant(L1,L2) - (b*c-a)^2*(b-c)^2/e2^2;
printf "Res(L2,L3) - (b*c-a)^2 = %o\n", Resultant(L2,L3) - (b*c-a)^2;
printf "L1(r)  - (b*c-a)*(b0-t)^2/e2 = %o\n", Evaluate(L1, r)  - (b*c-a)*(b0-t)^2/e2;
printf "L1(-r) - (b*c-a)*(b0+t)^2/e2 = %o\n", Evaluate(L1,-r) - (b*c-a)*(b0+t)^2/e2;
printf "L3(r) - (b*c-a) = %o\n", Evaluate(L3,r) - (b*c-a);

printf "=== (3) c-b == Psi modulo squares ===\n";
Psi := v*(1 + s^2*v)*(1 + A*v);
rat := (c-b)/Psi;
nm := Numerator(rat); dn := Denominator(rat);
printf "(c-b)/Psi = %o / %o\n", nm, dn;
okn, sq1 := IsSquare(nm*dn);
printf "IsSquare( num*den of (c-b)/Psi ) = %o\n", okn;
if okn then printf "   sqrt = %o\n", sq1; end if;

printf "=== (4) the halving system above T2 ===\n";
// u(x) = x^2 + w*x - (w*t^2 + b*c) is forced by u(r)*u(-r) being the Res(L2,.) value;
// verify the two 'value at +-r' conditions of  k u^2 + d2 L1 L3 = 0 mod L2.
Fw<w> := FunctionField(F);
Pw<X> := PolynomialRing(Fw);
L1w := Pw ! [ Fw!cc : cc in Coefficients(L1) ];
L2w := Pw ! [ Fw!cc : cc in Coefficients(L2) ];
L3w := Pw ! [ Fw!cc : cc in Coefficients(L3) ];
uw  := X^2 + w*X - (w*(Fw!t)^2 + Fw!(b*c));
kw  := -(Fw!d2)*(Fw!(b*c-a))^2/((Fw!e2)*(Fw!t)^2*w^2);
Rw  := kw*uw^2 + (Fw!d2)*L1w*L3w;
rem := Rw mod L2w;
printf "remainder of (k u^2 + d2 L1 L3) mod L2  =  %o\n", rem;
printf "lc(quotient) - (k+d2) = %o\n", LeadingCoefficient(Rw div L2w) - (kw + Fw!d2);

// the residual condition: the quotient must be a perfect square (A+Bx)^2
qq := Rw div L2w;
dsc := Coefficient(qq,1)^2 - 4*Coefficient(qq,0)*Coefficient(qq,2);
dn2 := Numerator(dsc);
printf "disc(quotient) numerator factorization:\n";
for g in Factorization(dn2) do
  printf "   [deg %o, mult %o] %o\n", Degree(g[1]), g[2],
         (#Sprint(g[1]) lt 900 select Sprint(g[1]) else "(large)");
end for;

printf "SYMBOLIC_DONE\n";
quit;
