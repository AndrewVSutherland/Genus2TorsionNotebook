//////////////////////////////////////////////////////////////////////
//  One-parameter M(12) family with simple examples and Z/12 x Z/2.
//
//  The extra rational Weierstrass point condition on M(12) has a simple
//  component:
//      a = (1-r)/4.
//  Then the completed-square polynomial factors as three rational
//  linear factors times a cubic, so the Jacobian contains an independent
//  rational 2-torsion point in addition to the order-12 point.
//
//  This script verifies the sample r=-25 and prints the family factors.
//////////////////////////////////////////////////////////////////////

Q := Rationals();
K<r> := RationalFunctionField(Q);
P<x> := PolynomialRing(K);

a := (1-r)/4;
T := a*x^2 - x + r;
h := (x-r)*(T+1);
f_rhs := a*x^2*T*(T+1);
W := h^2 + 4*f_rhs;
Q4 := W div (T+1);

print "a =", a;
print "T+1 factors:";
print Factorization(T+1);
print "W factors:";
print Factorization(W);
print "Q4 factors:";
print Factorization(Q4);

// Verify a concrete simple specialization.
Qx<X> := PolynomialRing(Q);
r0 := Q!-25;
a0 := Q!Evaluate(a, r0);
T0 := a0*X^2 - X + r0;
h0 := (X-r0)*(T0+1);
f0_rhs := a0*X^2*T0*(T0+1);
W0 := h0^2 + 4*f0_rhs;

assert Factorization(T0+1)[1][1] eq X - 2;

// Move w=2 to infinity.
w := Q!2;
f5 := Qx!0;
for i in [0..Degree(W0)] do
    for j in [0..i] do
        f5 +:= Coefficient(W0, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
    end for;
end for;
assert Degree(f5) eq 5;

C := HyperellipticCurve(f5);
J := Jacobian(C);
Xp := Q!-1/2;
Yp := Evaluate(h0,0)*Xp^3;
D := J![X-Xp, Yp];
assert Order(D) eq 12;

roots := Roots(f5);
assert #roots ge 2;
independent_found := false;
for rt in roots do
    beta := rt[1];
    Tbeta := J![X-beta, Q!0];
    if Tbeta ne J!0 and Tbeta ne 6*D then
        independent_found := true;
        print "independent rational 2-torsion root:", beta;
    end if;
end for;
assert independent_found;

Lp := LPolynomial(ChangeRing(C, GF(11)));
print "sample r =", r0;
print "odd quintic f5 =", f5;
print "LPolynomial at 11 =", Lp;
assert #Factorization(Lp) eq 1 and Degree(Factorization(Lp)[1][1]) eq 4;
print "Verified sample: order 12, independent 2-torsion, Q-simple certificate at p=11.";
quit;
