//Written by GPT 5.4 Thinking (small edits in Section 3)
//This is a companion file to half_divisor_gpt54.tex
//////////////////////////////////////////////////////////////////////
//  1. Generic check on the transformed curve
//
//  Curve:  Y^2 = z(z+1)(z+rho^2)(z+sigma^2)(z+tau^2)
//  Claim:  [ z^2 - s2*z + s4 , (s1*s2 - s3)*z - s1*s4 ]
//          doubles to [z,0], i.e. to (0,0)-infinity.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
K<rho,sigma,tau> := RationalFunctionField(Q, 3);
Pz<z> := PolynomialRing(K);

s1 := 1 + rho + sigma + tau;
s2 := rho + sigma + tau + rho*sigma + rho*tau + sigma*tau;
s3 := rho*sigma + rho*tau + sigma*tau + rho*sigma*tau;
s4 := rho*sigma*tau;

ftr := z*(z+1)*(z+rho^2)*(z+sigma^2)*(z+tau^2);
Ctr := HyperellipticCurve(ftr);
Jtr := Jacobian(Ctr);

Utr := z^2 - s2*z + s4;
Vtr := (s1*s2 - s3)*z - s1*s4;

H := Jtr![Utr, Vtr];
T0 := Jtr![z, 0];

assert 2*H eq T0;
print "Section 1 passed: generic halving formula on transformed curve verified.";


//////////////////////////////////////////////////////////////////////
//  2. Generic symbolic check of the pullback to the original x-model
//
//  We verify:
//
//    (a) the birational change of variables
//    (b) the pulled-back U(x)
//    (c) the explicit expanded V(x)
//
//  This part does NOT use Jacobian arithmetic; it is just algebra.
//////////////////////////////////////////////////////////////////////

K2<A,B,C,D,bv,s1g,s2g,s3g,s4g> := RationalFunctionField(Q, 9);

// ----- (a) Check the curve transformation algebraically -----

L<zz,YY> := RationalFunctionField(K2, 2);

xx := -(A + B*zz)/(1 + zz);
yy := bv*(B - A)/(1 + zz)^3 * YY;

// cleared numerator of y^2 - x(x+A)(x+B)(x+C)(x+D)
num := Numerator(yy^2 - xx*(xx+A)*(xx+B)*(xx+C)*(xx+D));

// expected cleared identity
target :=
-(A - B)^2 *
(
    bv^2*YY^2
    - zz*(1+zz)*(A + B*zz)*(A - C + (B - C)*zz)*(A - D + (B - D)*zz)
);

assert num eq target;
print "Section 2(a) passed: birational change of variables checked.";


// ----- (b) Check the pulled-back U(x) -----

Px<x> := PolynomialRing(K2);

delta := 1 + s2g + s4g;

U_from_pullback :=
((x + A)^2 + s2g*(x + A)*(x + B) + s4g*(x + B)^2) / delta;

U_explicit :=
x^2
+ ((2*A + s2g*(A+B) + 2*B*s4g)/delta)*x
+ (A^2 + A*B*s2g + B^2*s4g)/delta;

assert U_from_pullback eq U_explicit;
print "Section 2(b) passed: explicit U(x) checked.";


// ----- (c) Check the expanded V(x) -----

m := s1g*s2g - s3g;

// pulled-back cubic before reduction mod U
Vpull :=
(bv/(B-A)^2) * (x+B)^2 * ( -m*(x+A) - s1g*s4g*(x+B) );

// divide by U_explicit and take remainder
q, Vrem := Quotrem(Vpull, U_explicit);

Lambda :=
    s1g*s2g - s1g*s4g^2 + 3*s1g*s4g + s2g*s3g*s4g + 3*s3g*s4g - s3g;

M :=
    A*s1g*s2g + 2*A*s1g*s4g + A*s3g*s4g - A*s3g
    - B*s1g*s4g^2 + B*s1g*s4g + B*s2g*s3g*s4g + 2*B*s3g*s4g;

Vexplicit := -(bv/delta^2) * (Lambda*x + M);

assert Vrem eq Vexplicit;
print "Section 2(c) passed: explicit V(x) checked.";


//////////////////////////////////////////////////////////////////////
//  3. Optional concrete sanity check on the original curve
//
//  Example:
//      a = 17, b = 28, c = 32, d = 68
//
//  Then
//      u^2 = (a^2-c^2)(a^2-d^2) = 1785^2
//      v^2 = (b^2-c^2)(b^2-d^2) =  960^2
//      w^2 = (a^2-c^2)(b^2-c^2) =  420^2
//      t^2 = (a^2-d^2)(b^2-d^2) = 4080^2
//
//  We check directly in the Jacobian of the original curve that
//      2 D' = (P_a - inf) + (P_b - inf).
//////////////////////////////////////////////////////////////////////

Qx := Rationals();

//TO DO: replace with a loop through all of the 4-tuples in tor2244.txt, where each list of 4 elements is [a,b,c,d]
a := Qx!17;
b := Qx!28;
c := Qx!32;
d := Qx!68;

A0 := a^2;  B0 := b^2;  C0 := c^2;  D0 := d^2;
//u0 := Sqrt((A0-C0)*(A0-D0));
//v0 := Sqrt((B0-C0)*(B0-D0));
//w0 := Sqrt((A0-C0)*(B0-C0));
//t0 := Sqrt((A0-D0)*(B0-D0));
u0 := Qx!1785;
v0 := Qx!960;
w0 := Qx!420;
t0 := Qx!4080;

assert u0^2 eq (A0-C0)*(A0-D0);
assert v0^2 eq (B0-C0)*(B0-D0);
assert w0^2 eq (A0-C0)*(B0-C0);
assert t0^2 eq (A0-D0)*(B0-D0);
assert u0*v0 eq w0*t0;

rho0 := a/b;
sigma0 := (A0-C0)/w0;   // = w0/(B0-C0)
tau0 := (A0-D0)/t0;     // = t0/(B0-D0)

s10 := 1 + rho0 + sigma0 + tau0;
s20 := rho0 + sigma0 + tau0 + rho0*sigma0 + rho0*tau0 + sigma0*tau0;
s30 := rho0*sigma0 + rho0*tau0 + sigma0*tau0 + rho0*sigma0*tau0;
s40 := rho0*sigma0*tau0;
delta0 := 1 + s20 + s40;

Px0<x0> := PolynomialRing(Qx);
f0 := x0*(x0+A0)*(x0+B0)*(x0+C0)*(x0+D0);
C0curve := HyperellipticCurve(f0);
J0 := Jacobian(C0curve);

U0 :=
x0^2
+ ((2*A0 + s20*(A0+B0) + 2*B0*s40)/delta0)*x0
+ (A0^2 + A0*B0*s20 + B0^2*s40)/delta0;

Lambda0 :=
    s10*s20 - s10*s40^2 + 3*s10*s40 + s20*s30*s40 + 3*s30*s40 - s30;

M0 :=
    A0*s10*s20 + 2*A0*s10*s40 + A0*s30*s40 - A0*s30
    - B0*s10*s40^2 + B0*s10*s40 + B0*s20*s30*s40 + 2*B0*s30*s40;

V0 := -(b*v0/delta0^2) * (Lambda0*x0 + M0);

Dprime := J0![U0, V0];

// target class = (P_a - inf) + (P_b - inf)
Dtarget := J0![x0 + A0, 0] + J0![x0 + B0, 0];

assert 2*Dprime eq Dtarget;
print "Section 3 passed: concrete check on the original curve verified.";