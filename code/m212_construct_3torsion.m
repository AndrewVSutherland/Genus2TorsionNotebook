//////////////////////////////////////////////////////////////////////
//  Algebraic 3-torsion on the M(2,12) family.
//
//  Start from the M(12) model
//      y^2 + (x-r)(T+1)y = a*x^2*T*(T+1),  T=a*x^2-x+r.
//
//  On the M(2,12) chart, T+1 has a rational root.  Write
//      1 - 4*a*(r+1) = z^2,
//      a = (1-z^2)/(4*(r+1)).
//
//  We move the root w=2*(r+1)/(1+z) of T+1 to infinity.  The image of
//  P=(0,0) gives a divisor class D of order 12, so 4D is a rational
//  3-torsion class.  This script prints a generic Mumford representative.
//
//  Typical run:
//      magma code/m212_construct_3torsion.m
//////////////////////////////////////////////////////////////////////

SetColumns(0);

Q := Rationals();
K<z,r> := RationalFunctionField(Q, 2);
P<X> := PolynomialRing(K);

function OddQuinticAtRoot(W, w)
    out := P!0;
    for i in [0..Degree(W)] do
        for j in [0..i] do
            out +:= Coefficient(W, i)*Binomial(i,j)*w^(i-j)*X^(6-j);
        end for;
    end for;
    return P!out;
end function;

a := (1-z^2)/(4*(r+1));
T := a*X^2 - X + r;
h := (X-r)*(T+1);
W := h^2 + 4*a*X^2*T*(T+1);
w := 2*(r+1)/(1+z);

assert Evaluate(T+1, w) eq 0;

f5 := OddQuinticAtRoot(W, w);
Xp := -1/w;
Yp := Evaluate(h, 0)*Xp^3;

print "M(2,12) generic odd model after moving w to infinity";
print "a =", a;
print "w =", w;
print "f5 =", f5;
print "P_image =", <Xp, Yp>;

C := HyperellipticCurve(f5);
J := Jacobian(C);
D := J![X-Xp, Yp];
D3 := 4*D;

print "D =", D;
print "4D =", D3;
print "3*(4D) =", 3*D3;
print "nonzero_check_symbolic", D3 ne J!0;

quit;
